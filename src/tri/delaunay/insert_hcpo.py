'''
Created on Nov 13, 2018

@author: martijn
'''
import logging
import time
from math import ceil, sqrt
from random import shuffle, randint, random
from datetime import datetime
from operator import itemgetter
from tri.delaunay.tds import box, ccw, apex, orig, dest, \
    InfiniteVertex, Triangle, Vertex
from tri.delaunay.preds import orient2d, incircle
from tri.delaunay.tds import Triangulation
from tri.delaunay.iter import FiniteEdgeIterator
from tri.delaunay.cdt import ConstraintInserter
from tri.delaunay.inout import output_triangles, output_vertices
# py3k compatibility fix
try:
    # py2
    xrange
except NameError:
    # py3
    xrange = range


# -----------------------------------------------------------------------------
# Delaunay triangulation using Lawson's incremental insertion
#
def triangulate(points, infos=None, segments=None, output=False):
    """Triangulate a list of points, and if given also segments are
    inserted in the triangulation.
    """
    # FIXME: also embed info for points, if given as 3rd value in tuple
    # for every point
    if len(points) == 0:
        raise ValueError("we cannot triangulate empty point list")
    logging.debug("start " + str(datetime.now()))
    logging.debug("")
    logging.debug("pre-processing")
    start = time.clock()
    # points without info
    points = [(pt[0], pt[1], key) for key, pt in enumerate(points)]
    # this randomizes the points and then sorts them for spatial coherence
    points = hcpo(points)
    # get the original position and the new position in the sorted list
    # to build a lookup table for segment indices
    if infos is not None or segments is not None:
        index_translation = dict(((pos, newpos) for (newpos, (_, _, pos)) in
                                  enumerate(points)))
        if segments is not None:
            # -- translate the segments
            segments = [(index_translation[segment[0]],
                         index_translation[segment[1]])
                        for segment in segments]
        if infos is not None:
            infos = [(index_translation[idx], info)
                     for idx, info in enumerate(infos)]
    end = time.clock()
    logging.debug(str(end - start) + " secs")
    logging.debug("")
    logging.debug("triangulating " + str(len(points)) + " points")
    # add points, using incremental construction triangulation builder
    dt = Triangulation()
    start = time.clock()
    incremental = PointInserter(dt)
    incremental.insert(points)
    end = time.clock()
    logging.debug(str(end - start) + " secs")
    logging.debug(str(len(dt.vertices)) + " vertices")
    logging.debug(str(len(dt.triangles)) + " triangles")
    logging.debug(str(incremental.flips) + " flips")
    logging.debug(str(incremental.visits) + " visits")

    if len(dt.vertices) > 0:
        logging.debug(str(float(incremental.flips) /
                          len(dt.vertices)) + " flips per insert")

    # check links of triangles
#     check_consistency(dt.triangles)

    # insert segments
    if segments is not None:
        start = time.clock()
        logging.debug("")
        logging.debug("inserting " + str(len(segments)) + " constraints")
        constraints = ConstraintInserter(dt)
        constraints.insert(segments)
        end = time.clock()
        logging.debug(str(end - start) + " secs")
        logging.debug(str(len(dt.vertices)) + " vertices")
        logging.debug(str(len(dt.triangles)) + " triangles")
        edge_it = FiniteEdgeIterator(dt, constraints_only=True)
        constraint_ct = sum(1 for _ in edge_it)
        logging.debug(" {count} constraints".format(count=constraint_ct))
    # insert information for vertices
    if infos is not None:
        logging.debug("")
        logging.debug("inserting " + str(len(infos)) + " info")
        for info in infos:
            dt.vertices[info[0]].info = info[1]
    logging.debug("")
    logging.debug("fin " + str(datetime.now()))
    if output:
        with open("/tmp/alltris.wkt", "w") as fh:
            output_triangles(dt.triangles, fh)
        with open("/tmp/allvertices.wkt", "w") as fh:
            output_vertices(dt.vertices, fh)
        # FIXME: OUTPUT EDGES
    return dt


def cpo(points, c=0.5):
    """Column prime order for a set of points"""
    result = []
    (xmin, ymin), (xmax, ymax) = box(points)
    # width, height
    W = float(xmax - xmin)
    H = float(ymax - ymin)
    if H == 0:
        # prevents division by zero when calculating m
        H = 1.
    # sort on widest axis
    if W < H:
        axis = 1
    else:
        axis = 0
    points.sort(key=itemgetter(axis))
    # number of points to sort
    n = len(points)
    # determine bin size: how many bins do we need?
    m = int(ceil(c * ceil(sqrt(n * W / H))))
    if m == 0:
        # pathological case, no sampled points, so make it same as
        # number of points left
        m = n
    M = int(ceil(float(n) / float(m)))
    for i in range(m):
        j = i + 1
        # get bounds for this slot
        f, t = i * M, j * M
        slot = points[f:min(t, n)]
        # sort on other axis, in case even slot, in reversed order
        even = (j % 2) == 0
        slot.sort(key=itemgetter((axis + 1) % 2), reverse=even)

        # taking the same axis for sorting
#         slot.sort(key=itemgetter(axis), reverse=even) # twice as slow
        result.extend(slot)
    return result


def _hcpo(points, out, sr=0.75, minsz=10):
    """Constructs a hierarchical set of ordered points `out'

    Every level in the hierarchy is ordered along a column prime curve,
    similar to:

    >-------------------+
                        |
    +-------------------+
    |
    +-------------------+
                        |
    <-------------------+
    """
    # shuffle(points) # always randomize even for small points
    stack = [points]
    while stack:
        # split the remaining list in 2 pieces
        # tail will be processed (sorted and cut into slots)
        # head will be recursed on if it has enough points
        points = stack.pop()
        N = len(points)
        up = int(ceil(N*sr))
        head, tail = points[0:up], points[up:]
        if tail:
            ordered = cpo(tail)
            out.extend(ordered)
        if len(head) >= ceil(minsz / sr):
            shuffle(head)
            stack.append(head)
        else:
            ordered = cpo(head)
            out.extend(ordered)


def hcpo(points, sr=0.75, minsz=10):
    """Based on list with points, return a new, randomized list with points
    where the points are randomly ordered, but then sorted with enough spatial
    coherence to be useful to not get worst case flipping behaviour
    """
    # Build a new list with points, ordered along hierarchical curve
    # with levels
    if len(points) == 0:
        raise ValueError("not enough points")
    out = []
    _hcpo(points, out, sr, minsz)
    return out


class PointInserter(object):
    """Class to insert points into a Triangulation.

    It is ensured that the triangles that are made, are obeying the Delaunay
    criterion by flipping (Lawson's incremental algorithm is used
    to construct the triangulation).
    """

    __slots__ = ('triangulation', 'queue', 'flips', 'visits', 'last')

    def __init__(self, triangulation):
        self.triangulation = triangulation
        self.flips = 0
        self.visits = 0
        self.queue = []
        self.last = None  # last triangle used for finding triangle

    def insert(self, points):
        """Insert a list of points into the triangulation.
        """
        self.initialize(points)
        for j, pt in enumerate(points):
            self.append(pt)
            if (j % 10000) == 0:
                logging.debug(" " + str(datetime.now()) + str(j))
            # check_consistency(triangles)

    def initialize(self, points):
        """Initialize large triangle around point and external / dummy triangle
        from where we can always start point location
        """
        (xmin, ymin), (xmax, ymax) = box(points)
        width = abs(xmax - xmin)
        height = abs(ymax - ymin)
        if height > width:
            width = height
        if width == 0:
            width = 1.
        vertices = [InfiniteVertex(xmin - 50.0 * width, ymin - 40.0 * width),
                    InfiniteVertex(xmax + 50.0 * width, ymin - 40.0 * width),
                    InfiniteVertex(0.5 * (xmin + xmax), ymax + 60.0 * width)]
        large = Triangle(vertices[0], vertices[1], vertices[2])
        triangles = self.triangulation.triangles
        triangles.append(large)
        for v in vertices:
            v.triangle = large

    def append(self, pt):
        """Appends one point to the triangulation.

        This method assumes that the triangulation is initialized
        and the point lies inside the bounding box used for initializing.
        """
        v = Vertex(pt[0], pt[1])
        t0 = self.get_triangle_contains(v)
        # skip insertion of point, if it is on same location already there
        for corner in t0.vertices:
            if corner.x == v.x and corner.y == v.y:
                raise ValueError("Duplicate point found for insertion")
        self.triangulation.vertices.append(v)
        a, b, c = t0.vertices
        # neighbours outside triangle to insert to
        neighbours = [t0.neighbours[0], t0.neighbours[1]]
        neighbouridx = [n.neighbours.index(
            t0) if n is not None else None for n in neighbours]
        # make new triangles
        t1 = Triangle(b, c, v)
        t2 = Triangle(c, a, v)
        t0.vertices[2] = v
        # update triangle pointers of vertices
        a.triangle = t0
        b.triangle = t0
        v.triangle = t0
        c.triangle = t1
        # link them up properly -- use neighbours outside triangle to insert to
        # external links
        # 2 * 2
        if neighbours[0] is not None:
            side = neighbouridx[0]
            self.link_1dir(neighbours[0], side, t1)
        self.link_1dir(t1, 2, neighbours[0])
        if neighbours[1] is not None:
            side = neighbouridx[1]
            self.link_1dir(neighbours[1], side, t2)
        self.link_1dir(t2, 2, neighbours[1])
        # internal links
        # 3 * 2
        self.link_2dir(t0, 0, t1, 1)
        self.link_2dir(t1, 0, t2, 1)
        self.link_2dir(t2, 0, t0, 1)
        #
        triangles = self.triangulation.triangles
        triangles.extend([t1, t2])
        # check if triangles are delaunay, and flip
        # edges of triangle just inserted into are queued for checking
        # Delaunay criterion
        self.queue.append((t2, 2))
        self.queue.append((t1, 2))
        self.queue.append((t0, 2))
        self.delaunay()

    def get_triangle_contains(self, p):
        """Gets the triangle on which point p is located from the triangulation
        """
        # ini = self.random_triangle_close_to_p(p)
        # FIXME:
        if self.triangulation.vertices:
            ini = self.triangulation.vertices[-1].triangle
        else:
            ini = self.random_triangle_close_to_p(p)
        # Pick random close triangle or
        # start walk from last inserted point triangle
        t0 = self.visibility_walk(ini, p)
        # remember this triangle as it might be close to next wanted point
        self.last = t0
        return t0

    def random_triangle_close_to_p(self, p):
        """Samples a list of triangles and returns closest of these triangles
        to the given point p
        """
        # FIXME: should we cache result of random triangles
        # as long as sample size is the same
        # we could use the same set of triangles
        # O(n/3) would be good where n is the number of triangles
        candidate = self.triangulation.triangles[0]  # external
        min_dist = None  # candidate.vertices[0].distance(p)
        triangles = self.triangulation.triangles
        size = len(triangles)
        #
        if size != 0:
            k = int(sqrt(size) / 25)
            # k = int(size ** (1 / 3.0)) # -- samples more triangles
            if self.last is not None:  # set by triangle_contains
                dist = self.last.vertices[0].distance2(p)
                if min_dist is None or dist < min_dist:
                    min_dist = dist
                    candidate = self.last
            for _ in range(k):
                triangle = triangles[int(random() * size)]
                dist = triangle.vertices[0].distance2(p)
                if min_dist is None or dist < min_dist:
                    min_dist = dist
                    candidate = triangle
        return candidate

    def visibility_walk(self, ini, p):
        """Walk from triangle ini to triangle containing p

        Note, because this walk can cycle for a non-Delaunay triangulation
        we pick a random edge to continue the walk
        (this is a remembering stochastic walk, see RR-4120.pdf,
        Technical report from HAL-Inria by
        Olivier Devillers, Sylvain Pion, Monique Teillaud.
        Walking in a triangulation,
        https://hal.inria.fr/inria-00072509)

        For speed we do not check if we stay inside the bounding box
        that was used when initializing the triangulation, so make sure
        that a point given fits inside this box!
        """
        t = ini
        previous = None
        if t.vertices[2] is None:
            t = t.neighbours[2]
        n = len(self.triangulation.triangles)
        for ct in xrange(n):
            # get random side to continue walk, this way the walk cannot get
            # stuck by always picking triangles in the same order
            # (and get stuck in a cycle in case of non-Delaunay triangulation)
            e = randint(0, 2)
            if t.neighbours[e] is not previous and \
                orient2d(t.vertices[ccw(e)],
                         t.vertices[ccw(e+1)],
                         p) < 0:
                previous = t
                t = t.neighbours[e]
                continue
            e = ccw(e + 1)
            if t.neighbours[e] is not previous and \
                orient2d(t.vertices[ccw(e)],
                         t.vertices[ccw(e+1)],
                         p) < 0:
                previous = t
                t = t.neighbours[e]
                continue
            e = ccw(e + 1)
            if t.neighbours[e] is not previous and \
                orient2d(t.vertices[ccw(e)],
                         t.vertices[ccw(e+1)],
                         p) < 0:
                previous = t
                t = t.neighbours[e]
                continue
            self.visits += ct
            return t
        self.visits += ct
        return t

    def delaunay(self):
        """Flips triangles if Delaunay criterion does not hold.

        If 2 triangles were flipped, the 4 triangles around the quadrilateral
        are queued for checking if these are Delaunay.
        """
        while self.queue:
            t0, side0 = self.queue.pop()
            # -- skip constrained edge - these should not be flipped
            if t0.constrained[side0]:
                continue
            t1 = t0.neighbours[side0]
            # -- skip if we are going to flip the external dummy triangle
            # or when the triangle is an infinite triangle
#             if t1 is self.triangulation.external or t1 is None:
            if t1 is None:
                continue
            # -- get the opposite vertex/side index
            # it's an error if we cannot find t0
            side1 = t1.neighbours.index(t0)
            if side1 is None:
                raise ValueError("No opposite triangle found")
            if incircle(t0.vertices[0], t0.vertices[1], t0.vertices[2],
                        t1.vertices[side1]) > 0:
                # flip triangles without creating new triangle objects
                self.flip(t0, side0, t1, side1)
                # check if all 4 edges around quadrilateral just flipped
                # are now good: i.e. delaunay criterion applies
                self.queue.append((t0, 0))
                self.queue.append((t0, 2))
                self.queue.append((t1, 0))
                self.queue.append((t1, 2))

    def flip(self, t0, side0, t1, side1):
        """Performs the flip of triangle t0 and t1

        If t0 and t1 are two triangles sharing a common edge AB,
        the method replaces ABC and BAD triangles by DCA and DBC, respectively.

        Pre-condition:
        triangles t0/t1 share a common edge and the edge is known
        """
        self.flips += 1

        apex0, orig0, dest0 = apex(side0), orig(side0), dest(side0)
        apex1, orig1, dest1 = apex(side1), orig(side1), dest(side1)

        # side0 and side1 should be same edge
        assert t0.vertices[orig0] is t1.vertices[dest1]
        assert t0.vertices[dest0] is t1.vertices[orig1]

        # assert both triangles have this edge unconstrained
        assert not t0.constrained[apex0]
        assert not t1.constrained[apex1]

        # -- vertices around quadrilateral in ccw order starting at apex of t0
        A, B = t0.vertices[apex0], t0.vertices[orig0]
        C, D = t1.vertices[apex1], t0.vertices[dest0]
        # -- triangles around quadrilateral in ccw order, starting at A
        AB, BC = t0.neighbours[dest0], t1.neighbours[orig1]
        CD, DA = t1.neighbours[dest1], t0.neighbours[orig0]
        # link neighbours around quadrilateral to triangles as after the flip
        # -- the sides of the triangles around are stored in apex_around
        apex_around = []
        for neighbour, corner in zip([AB, BC, CD, DA],
                                     [A, B, C, D]):
            if neighbour is None:
                apex_around.append(None)
            else:
                apex_around.append(ccw(neighbour.vertices.index(corner)))
        # the triangles around we link to the correct triangle *after* the flip
        for neighbour, side, t in zip([AB, BC, CD, DA],
                                      apex_around,
                                      [t0, t0, t1, t1]):
            if neighbour is not None:
                self.link_1dir(neighbour, side, t)

        # -- set new vertices and neighbours
        # for t0
        t0.vertices = [A, B, C]
        t0.neighbours = [BC, t1, AB]
        # for t1
        t1.vertices = [C, D, A]
        t1.neighbours = [DA, t0, CD]
        # -- update coordinate to triangle pointers
        for v in t0.vertices:
            v.triangle = t0
        for v in t1.vertices:
            v.triangle = t1

    def link_2dir(self, t0, side0, t1, side1):
        """Links two triangles to each other over their common side
        """
        assert t0 is not None
        assert t1 is not None
        t0.neighbours[side0] = t1
        t1.neighbours[side1] = t0

    def link_1dir(self, t0, side0, t1):
        """Links triangle t0 to t1 for side0"""
        t0.neighbours[side0] = t1


def check_consistency(triangles):
    """Check a list of triangles for consistent neighbouring relationships

    For every triangle in the list
    it checks whether the triangle its neighbours also
    point back to this triangle.
    """
    errors = []
    for t in triangles:
        for n in t.neighbours:
            if n is not None:
                if t not in n.neighbours:
                    errors.append("{} {}".format(id(t), id(n)))
    if len(errors) > 0:
        raise ValueError("\n".join(errors))
