'''
Created on Nov 14, 2018

@author: martijn
'''

import operator
import logging
import time
from datetime import datetime
from random import randint

from tri.delaunay.tds import box, ccw, apex, orig, dest, InfiniteVertex, \
    Vertex, Triangle, Triangulation
from tri.delaunay.preds import orient2d, incircle
from tri.delaunay.cdt import ConstraintInserter
from tri.delaunay.iter import FiniteEdgeIterator, StarEdgeIterator
from tri.delaunay.inout import output_vertices, output_triangles

# py3k
try:
    xrange
except NameError:
    xrange = range


def decorate(points):
    """Adds index to every item in points list
    (every item is dealt with as 2-tuple)

    Returns a list with 3-tuples, where every 3-tuple contains:

        (x, y, index in the original *points* list)
    """
    ret = [(pt[0], pt[1], idx) for (idx, pt) in enumerate(points, start=0)]
    return ret


def largest_axis(aabb):
    """Given an axis-aligned bounding box as two 2-tuples, what is the
    largest axis of this box
    """
    dx = aabb[1][0] - aabb[0][0]
    dy = aabb[1][1] - aabb[0][1]
    if dx > dy:
        return 0
    else:
        return 1


def kdsort(points):
    """Sorts a list of tuples based on first two elements along kD-tree order.

    To every tuple in the result list the parent point index (according
    to the kD-tree order) is added, so that this point can be used
    to start a walk in the triangulation from that vertex.
    """
    stack = []  # low, high, parent, box
    result = []
    stack.append((points[:], None, box(points)))
    while stack:
        points, parent_id, aabb = stack.pop()
        axis = largest_axis(aabb)
        points.sort(key=operator.itemgetter(axis))
        halfway = len(points) // 2
        # get the pivot point and add the parent vertex identifier to it
        pivot = tuple(list(points[halfway]) + [parent_id])
        result.append(pivot)
        # determine the next halves
        (xmid, ymid) = pivot[0], pivot[1]
        (left, bottom) = aabb[0]
        (right, top) = aabb[1]
        leftpts = points[:halfway]
        rightpts = points[halfway+1:]
        # stack right half
        if rightpts:
            if axis == 0:
                half_aabb = [(xmid, bottom), (right, top)]
            elif axis == 1:
                half_aabb = [(left, ymid), (right, top)]
            stack.append(
                (rightpts, len(result) - 1, half_aabb)
            )
        # stack left half
        if leftpts:
            if axis == 0:
                half_aabb = [(left, bottom), (xmid, top)]
            elif axis == 1:
                half_aabb = [(left, bottom), (right, ymid)]
            stack.append(
                (leftpts, len(result) - 1, half_aabb)
            )
    return result


def translate_old2new(old_decorated_sorted_pts):
    """Build a translation table.
    """
    return dict([(pt[2], new_pos) for (new_pos, pt) in enumerate(old_decorated_sorted_pts)])

def translate_new2old(old_decorated_sorted_pts):
    """Build a translation table.
    """
    return dict([(new_pos, pt[2]) for (new_pos, pt) in enumerate(old_decorated_sorted_pts)])


class KDOrderPointInserter(object):
    """Class to insert points into a Triangulation.

    It is ensured that the triangles that are made, are obeying the Delaunay
    criterion by flipping (Lawson's incremental algorithm is used
    to construct the triangulation).
    """

    __slots__ = ('triangulation', 'queue', 'flips', 'visits', '_walk')

    def __init__(self, triangulation):
        self.triangulation = triangulation
        self.flips = 0
        self.visits = 0
        self.queue = []

        self._walk = self.visibility_walk

    def insert(self, points):
        """Insert a list of points into the triangulation.
        """
        self.initialize(points)
        for j, pt in enumerate(points):
            logging.debug(" - inserting {}".format(pt))
            parent = pt[3]
            ini = self.triangulation.triangles[0]
            if parent is not None:
                ini = self.triangulation.vertices[parent].triangle
            assert ini is not None
            self.append(pt, ini)
            if (j % 10000) == 0:
                logging.debug(" " + str(datetime.now()) + str(j))
            # check_consistency(triangles)

#    def initialize(self, points):
#        """Initialize large triangle around point and external / dummy triangle
#        from where we can always start point location
#        """
#        (xmin, ymin), (xmax, ymax) = box(points)
#        width = abs(xmax - xmin)
#        height = abs(ymax - ymin)
#        if height > width:
#            width = height
#        if width == 0:
#            width = 1.
#        vertices = [InfiniteVertex(xmin - 50.0 * width, ymin - 40.0 * width),
#                    InfiniteVertex(xmax + 50.0 * width, ymin - 40.0 * width),
#                    InfiniteVertex(0.5 * (xmin + xmax), ymax + 60.0 * width)]
#        large = Triangle(vertices[0], vertices[1], vertices[2])
#        self.triangulation.external = Triangle(vertices[1], vertices[0], None)
#        triangles = self.triangulation.triangles
#        triangles.append(large)
#        self.link_2dir(large, 2, self.triangulation.external, 2)
#        for v in vertices:
#            v.triangle = large


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
        inf_vertices = \
            [InfiniteVertex(xmin - 50.0 * width, ymin - 40.0 * width),
             InfiniteVertex(xmax + 50.0 * width, ymin - 40.0 * width),
             InfiniteVertex(0.5 * (xmin + xmax), ymax + 60.0 * width),
             None
             # InfiniteVertex(xmin + 0.5 * width, ymin + 0.5 * width)
             #
             # FIXME:
             # put also geometric embedding of this infinite vertex
             # and make algorithms robust
             # to not walk into external triangle!
             #
             # InfiniteVertex(xmin + 0.5 * width, ymin + 0.5 * width)
             # lifted tip of cone (in center of the domain...)
             ]
        large = Triangle(inf_vertices[0], inf_vertices[1], inf_vertices[2])
        hat = [Triangle(inf_vertices[1], inf_vertices[0], inf_vertices[3]),
               Triangle(inf_vertices[2], inf_vertices[1], inf_vertices[3]),
               Triangle(inf_vertices[0], inf_vertices[2], inf_vertices[3])]

        hat0 = hat[0]
        hat1 = hat[1]
        hat2 = hat[2]

        self.link_2dir(large, 2, hat0, 2)
        self.link_2dir(large, 0, hat1, 2)
        self.link_2dir(large, 1, hat2, 2)

        self.link_2dir(hat0, 1, hat1, 0)
        self.link_2dir(hat1, 1, hat2, 0)
        self.link_2dir(hat2, 1, hat0, 0)


#        self.triangulation.external = Triangle(inf_vertices[1], inf_vertices[0], None)
# FIXME: Is this 'external' triangle needed???

        # FIXME: Make 3 triangles around large triangle?
#         e0 = Triangle(inf_vertices[1], inf_vertices[0], None)
#         self.link_2dir(large, 2, e0, 2)
#         e1 = Triangle(inf_vertices[2], inf_vertices[1], None)
#         self.link_2dir(...)
#         e2 = Triangle(inf_vertices[0], inf_vertices[2], None)
#         self.link_2dir(...)

        triangles = self.triangulation.triangles
        triangles.append(large)
#        triangles.extend(hat)

#        self.link_2dir(large, 2, self.triangulation.external, 2)
        for v in inf_vertices:  # exclude the tip of the cone
            if v:  # skip None
                v.triangle = large

    def append(self, pt, ini):
        """Appends one point to the triangulation.

        This method assumes that the triangulation is initialized
        and the point lies inside the bounding box used for initializing.
        """
        v = Vertex(pt[0], pt[1])
        t0 = self._walk(ini, v)
        # skip insertion of point, if it is on same location already there
        for corner in t0.vertices:
            if corner.x == v.x and corner.y == v.y:
                raise ValueError("Duplicate point found for insertion")
        self.triangulation.vertices.append(v)
        a, b, c = t0.vertices
        # neighbours outside triangle to insert to
        neighbours = [t0.neighbours[0], t0.neighbours[1]]
        neighbouridx = [n.neighbours.index(t0)
                        if n is not None else None
                        for n in neighbours]
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
        self.triangulation.triangles.extend([t1, t2])
        # check if triangles are delaunay, and flip
        # edges of triangle just inserted into are queued for checking
        # Delaunay criterion
        self.queue.append((t2, 2))
        self.queue.append((t1, 2))
        self.queue.append((t0, 2))
        self.delaunay()


    def remove(self, pt, ini):
        """Remove a point that is present in the triangulation,
        while keeping the triangulation Delaunay.

        The algorithm followed is from the paper by
        Mostafavi, Gold & Dakowicz (2003).

        @article{Mostafavi2003,
          doi = {10.1016/s0098-3004(03)00017-7},
          url = {https://doi.org/10.1016/s0098-3004(03)00017-7},
          year = {2003},
          month = may,
          publisher = {Elsevier {BV}},
          volume = {29},
          number = {4},
          pages = {523--530},
          author = {Mir Abolfazl Mostafavi and Christopher Gold and Maciej Dakowicz},
          title = {Delete and insert operations in Voronoi/Delaunay methods and applications},
          journal = {Computers {\&} Geosciences}
        }

        Note 1: the removal does happily remove points if a CDT is used (not ok).

        Note 2: infinite vertex removal should not happen

        """
        ### FIXME: This removal does not care about presence of constraints in the DT!
        ### Maybe we should make the Triangulation structure aware of this
        ### (either by subclassing or by having a bool attribute)
        tri = self._walk(ini, Vertex(pt[0], pt[1]))
        pivot = None
        for v in tri.vertices:
            if (v[0] == pt[0] and v[1] == pt[1]):
                pivot = v
                break
        del v
        if pivot is None:
            raise ValueError("{} not found in triangulation, are you sure it was inserted?".format(pt))

        # -- slice ears for the polygon around the pivot that is to be removed
        # the polygon around the pivot to be removed is represented by a
        # collection of triangles, which we call the *star*
        # the vertices on this polygon are called the *link*
        star = [edge.triangle for edge in StarEdgeIterator(pivot)]
        cur = 0
        while len(star) > 3:
             # take 2 triangles (going around ccw around the pivot)
            tri0 = star[(cur) % len(star)]
            tri1 = star[(cur + 1) % len(star)]

            # get the vertices opposite of the pivot
            # tri0
            side0 = tri0.vertices.index(pivot)
            v1 = tri0.vertices[orig(side0)]
            v2 = tri0.vertices[dest(side0)]
            # tri1
            side1 = tri1.vertices.index(pivot)
            v_ = tri1.vertices[orig(side1)]
            v3 = tri1.vertices[dest(side1)]
            # they should share 1 vertex
            assert v2 is v_

            #print(" ear to check --> {}".format( 
            #    " | ".join(map(lambda v: "{:8.2f} {:8.2f}".format(v.x, v.y),
            #                   [v1, v2, v3]))))
            tri2 = star[(cur + 2) % len(star)]

            # we have a potential ear to slice off
            # if (v1, v2, v3 turns left) and (v1, v3, pivot turns left or is straight) 
            if orient2d(v1, v2, v3) > 0 and orient2d(v1, v3, pivot) >= 0:
                # make sure we really can slice the ear off 
                slice_ear = True
                # check for other triangles (except tri0, tri1, tri2) its orig()-point
                # (i.e. circumcircle through ear its points is empty of other points in the link)
                for i in range(0, len(star) - 3):
                    tri3 = star[(cur + 3 + i) % len(star)]
                    assert tri3 is not tri0
                    assert tri3 is not tri1
                    assert tri3 is not tri2
                    v4 = tri3.vertices[ orig(tri3.vertices.index(pivot)) ]
                    if incircle(v1, v2, v3, v4) > 0: # v4 inside circle, do not slice
                        slice_ear = False
                        break

                if slice_ear:
                    # -- Ear will be sliced, 
                    # flipping 2 triangles reduces the star by one triangle
                    #
                    # flip22 flips CCW
                    # -> tri0 is thus the sliced ear,
                    #    so remove tri0 from the star and keep only tri1 in the remaining star
                    self.flip22(tri0, tri0.vertices.index(v1),
                                tri1, tri1.vertices.index(v3))
                    star.remove(tri0)
                    # pivot should not be in tri0 its vertices
                    assert pivot not in tri0.vertices
                    # yet it should appear in tri1
                    assert pivot in tri1.vertices
#                    if raw_input('#! write [y/N] $ ') in ('y', 'Y'):
#                        with open("/tmp/all_tris.wkt", "w") as fh:
#                            print(len(self.triangulation.triangles))
#                            output_triangles(self.triangulation.triangles, fh)
#                        with open("/tmp/all_vertices.wkt", "w") as fh:
#                            output_vertices(self.triangulation.vertices, fh)
            # increment
            cur += 1
            # and -- possibly -- wrap around to start of list (make cur = 0)
            cur %= len(star)

#        with open("/tmp/tri__all_tris__before_flip31.wkt", "w") as fh:
#            output_triangles(self.triangulation.triangles, fh)
#        with open("/tmp/tri__all_vertices__before_flip31.wkt", "w") as fh:
#            output_vertices(self.triangulation.vertices, fh)
#        print(id(pivot.triangle), "to be removed")
        assert len(star) == 3
        # -- now remove the 3 triangles by performing a flip3->1
        self.flip31(pivot)


    def flip31(self, vertex):
        """ 'Flips' 3 triangles into 1,
        i.e. dissolves the three triangles around the pivot *vertex* into 1 triangle

        Pre-condition: the vertex is surrounded by exactly 3 triangles,
        this condition is _not_ checked
        """
        # given a vertex
        # remove 2 of its adjacent triangles that also share an edge with this vertex
        # keeping 1 tri in the triangulation,
        # linking it to the two neighbours of 
        # the 2 adjacent triangles that are removed
        # and hereby reducing the total number of triangles by two
        #
        # Note, these 2 triangles are garbage in the triangles list of the dt
        # and need to be garbage collected... similar to the vertex removed

        tri = vertex.triangle
        side0 = tri.vertices.index(vertex)
        apex0, orig0, dest0 = apex(side0), orig(side0), dest(side0)

        # the two triangles that will be removed
        ngb_orig = tri.neighbours[orig0] # neighbour tri opposite orig point
        ngb_dest = tri.neighbours[dest0] # neighbour tri opposite dest point

        # the 2 'around' triangles -- to link tri with!
        link_orig = ngb_orig.neighbours[ngb_orig.vertices.index(vertex)]
        link_dest =  ngb_dest.neighbours[ngb_dest.vertices.index(vertex)]

        # the new 'apex' of the new triangle
        if link_orig:
            tip_orig = link_orig.vertices[dest(link_orig.vertices.index(tri.vertices[dest0]))]
            tri.vertices[apex0] = tip_orig
        if link_dest:
            tip_dest = link_dest.vertices[orig(link_dest.vertices.index(tri.vertices[orig0]))]
            tri.vertices[apex0] = tip_dest

        # extra check (in case we have both triangles present,
        # they should both refer to the same new apex)
        if link_orig and link_dest:
            assert tip_orig is tip_dest

        # link neighbours tri to link_orig
        if link_orig:
            self.link_2dir(tri, orig0, link_orig, orig(link_orig.vertices.index(tri.vertices[dest0])))
        else:
            self.link_1dir(tri, orig0, None)
        # link neighbours tri to link_dest
        if link_dest:
            self.link_2dir(tri, dest0, link_dest, dest(link_dest.vertices.index(tri.vertices[orig0])))
        else:
            self.link_1dir(tri, dest0, None)
        # fix triangle pointers of the vertices of the triangle
        for v in tri.vertices:
            v.triangle = tri 
        #
        ngb_orig.neighbours = None
        ngb_orig.vertices = None
        ngb_orig.constrained = None
        ngb_orig.info = None
        #
        ngb_dest.neighbours = None
        ngb_dest.vertices = None
        ngb_dest.constrained = None
        ngb_dest.info = None
        #
        vertex.x = None
        vertex.y = None
        vertex.info = None
        vertex.triangle = None
        # FIXME: this will be slow for larger triangulations :'-(
        # note we could also leave the removed triangle/vertex objects in the list
        # and garbage collect them later, either when
        # we get too much garbage, or when we run a function
        # that depends on the triangles / vertices list in the triangulation object
        # (e.g. the output_* funcs)
        # or we remove the storage of the triangles / vertices in the list
        # completely and leave it up to the garbage collector of python

        # hybrid: keep them in the list, remove them based on own function call,
        # e.g. clean()

        self.triangulation.triangles.remove(ngb_orig)
        self.triangulation.triangles.remove(ngb_dest)
        self.triangulation.vertices.remove(vertex)


#    def get_triangle_contains(self, p):
#        """Gets the triangle on which point p is located from the triangulation
#        """
#        # ini = self.random_triangle_close_to_p(p)
#        ## FIXME:
#        if self.triangulation.vertices:
#            ini = self.triangulation.vertices[-1].triangle
#        else:
#            ini = self.random_triangle_close_to_p(p)
#        ## Pick random close triangle or start walk from last inserted point
#        ## triangle
#        t0 = self.visibility_walk(ini, p)
#        # remember this triangle as it might be close to next wanted point
#        self.last = t0
#        return t0

#    def random_triangle_close_to_p(self, p):
#        """Samples a list of triangles and returns closest of these triangles
#        to the given point p
#        """
#        # FIXME: should we cache result of random triangles
#        # as long as sample size is the same
#        # we could use the same set of triangles
#        # O(n/3) would be good where n is the number of triangles
#        candidate = self.triangulation.external
#        min_dist = None # candidate.vertices[0].distance(p)
#        triangles = self.triangulation.triangles
#        size = len(triangles)
#        #
#        if size != 0:
#            k = int(sqrt(size) / 25)
#            #k = int(size ** (1 / 3.0)) # -- samples more triangles
#            if self.last is not None: # set by triangle_contains
#                dist = self.last.vertices[0].distance2(p)
#                if min_dist is None or dist < min_dist:
#                    min_dist = dist
#                    candidate = self.last
#            for _ in xrange(k):
#                triangle = triangles[int(random() * size)]
#                dist = triangle.vertices[0].distance2(p)
#                if min_dist is None or dist < min_dist:
#                    min_dist = dist
#                    candidate = triangle
#        return candidate

    def straight_walk(self, ini, p):
        """Walk straight from triangle ini to triangle containing p"""
        tri = ini
        q_idx, r_idx, l_idx = 0, 1, 2
        # initialize the walk by rotating ccw or cw
        q = tri.vertices[q_idx]
        if orient2d(tri.vertices[r_idx], tri.vertices[q_idx], p) < 0:
            while orient2d(tri.vertices[l_idx], tri.vertices[q_idx], p) < 0:
                neighbour = tri.neighbours[r_idx]
                q_idx = neighbour.vertices.index(q)
                r_idx = ccw(q_idx)
                l_idx = ccw(r_idx)
                tri = neighbour
                assert tri.vertices[q_idx] == q
        else:
            while True:
                neighbour = tri.neighbours[l_idx]
                q_idx = neighbour.vertices.index(q)
                r_idx = ccw(q_idx)
                l_idx = ccw(r_idx)
                tri = neighbour
                assert neighbour.vertices[q_idx] == q
                if orient2d(tri.vertices[r_idx], tri.vertices[q_idx], p) < 0:
                    break
        # perform the walk
        s_idx = q_idx
        while orient2d(p, tri.vertices[r_idx], tri.vertices[l_idx]) < 0.0:
            neighbour = tri.neighbours[s_idx]
            l_idx = neighbour.vertices.index(tri.vertices[l_idx])
            r_idx = ccw(l_idx)
            s_idx = ccw(r_idx)
            tri = neighbour
            # 'advance' 1 of the two sides of the neighbour triangle
            # by swapping either l or r with s
            if orient2d(tri.vertices[s_idx], q, p) < 0.0:
                s_idx, r_idx = r_idx, s_idx
            else:
                s_idx, l_idx = l_idx, s_idx
        return tri

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
        """Performs Flip22 for triangles if Delaunay criterion does not hold.

        If 2 triangles were flipped, the 4 triangles around the quadrilateral
        are queued for checking if these are Delaunay.
        """
        while self.queue:
            t0, side0 = self.queue.pop()
            # -- skip constrained edge - these should not be flipped
            if t0.constrained[side0]:
                continue
            t1 = t0.neighbours[side0]
            # -- skip if we are about to flip with a 'non-existent' triangle (outside of big triangle --> 'hat')
            if t1 is None or t1.vertices[2] is None:
                continue
            # -- get the opposite vertex/side index
            # it's an error if we cannot find t0
            side1 = t1.neighbours.index(t0)
            if side1 is None:
                raise ValueError("No opposite triangle found")
            if incircle(t0.vertices[0], t0.vertices[1], t0.vertices[2],
                        t1.vertices[side1]) > 0:
                # flip triangles without creating new triangle objects
                self.flip22(t0, side0, t1, side1)
                # check if all 4 edges around quadrilateral just flipped
                # are now good: i.e. delaunay criterion applies
                self.queue.append((t0, 0))
                self.queue.append((t0, 2))
                self.queue.append((t1, 0))
                self.queue.append((t1, 2))

    def flip22(self, t0, side0, t1, side1):
        """Performs the flip of triangle t0 and t1

        If t0 and t1 are two triangles sharing a common edge AB,
        the method replaces ABC and BAD triangles by DCA and DBC, respectively.

        Pre-condition:
        Triangles t0/t1 share a common edge and the edge is known

        Post-conditions:
        - t0 / t1 are rotated *ccw* 
        - t0 / t1 are linked correctly within the quad (vertices/neighbouring triangles) and wrt each other
        - the vertices point to the correct triangle
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
                assert neighbour is not None
                assert neighbour.vertices is not None
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


def triangulate(pts, infos=None, segments=None, output=False):
    """Triangulate a set of points
    """

    # FIXME: this function would be better as 4 sub-functions, that can be re-used, e.g.
    # _preprocess_vertices()
    #    pts = kdsort(decorate(pts))
    #    new2old_index = translation_table(pts)
    # _insert_vertices()
    # _insert_infos()
    # _insert_segments
    # _output_triangulation

    start = time.perf_counter()
    logging.debug("")
    logging.debug(list(enumerate(["{}".format(_)for _ in pts])))
#     orig_pts = pts[:]

    pts = kdsort(decorate(pts))
    old2new = translate_old2new(pts)
    new2old = translate_new2old(pts)
    end = time.perf_counter()
    logging.debug("Sorting points: " + str(end - start) + " secs")

    start = time.perf_counter()
    dt = Triangulation()
    incremental = KDOrderPointInserter(dt)
    incremental.insert(pts)
    end = time.perf_counter()

    logging.debug("Triangulating took: " + str(end - start) + " secs")
    logging.debug("{} triangles".format(len(dt.triangles)))
    logging.debug("{} vertices".format(len(dt.vertices)))
    logging.debug("{} flips".format(incremental.flips))
    logging.debug("{} visits".format(incremental.visits))
    if len(dt.vertices) > 0:
        logging.debug(str(float(incremental.flips) /
                          len(dt.vertices)) + " flips per insert")

    if output:
        with open("/tmp/all_tris.wkt", "w") as fh:
            output_triangles(dt.triangles, fh)
        with open("/tmp/all_vertices.wkt", "w") as fh:
            output_vertices(dt.vertices, fh)

    if infos:
        logging.debug("adding info on vertices")
        for new_idx, v in enumerate(dt.vertices):
            # translate index (after kD-sort)
            old_idx = new2old[new_idx]
            v.info = infos[old_idx]

    if output:
        with open("/tmp/all_tris.wkt", "w") as fh:
            output_triangles(dt.triangles, fh)
        with open("/tmp/all_vertices.wkt", "w") as fh:
            output_vertices(dt.vertices, fh)

    if segments:
        start = time.perf_counter()
        logging.debug("")
        logging.debug("inserting " + str(len(segments)) + " constraints")
        constraints = ConstraintInserter(dt)
        # translate indexes (after kD-sort) of segments to be inserted
        logging.debug(list(enumerate(["{}".format(_) for _ in dt.vertices])))
        logging.debug(segments)
#         orig_segments = segments[:]
#         for segment in orig_segments:
#             print("old", orig_pts[segment[0]], orig_pts[segment[1]])
#         print("")
        segments = [(old2new[segment[0]], old2new[segment[1]])
                    for segment in segments]
#         for segment in segments:
#             print("new", pts[segment[0]][:2], pts[segment[1]][:2])
        logging.debug(segments)
#         from random import shuffle
#         print(shuffle(segments))
        constraints.insert(segments)
        end = time.perf_counter()
        logging.debug(" {time} secs".format(time=(end-start)))
        logging.debug(" {vertex_count} vertices".format(
                        vertex_count=len(dt.vertices)))
        logging.debug(" {triangle_count} triangles".format(
                        triangle_count=len(dt.triangles)))
        # Keep FiniteEdgeIterator as iterator (do not read it to memory)
        edge_it = FiniteEdgeIterator(dt, constraints_only=True)
        constraint_ct = sum(1 for _ in edge_it)
        logging.debug(" {count} constraints".format(count=constraint_ct))

    if output:
        with open("/tmp/all_tris.wkt", "w") as fh:
            output_triangles(dt.triangles, fh)
        with open("/tmp/all_vertices.wkt", "w") as fh:
            output_vertices(dt.vertices, fh)

    return dt


def test_sort_keeps_connection_info():
    """The sorting process should maintain the connection to the information
    associated with each vertex.
    """
    pts = [(0., 0.), (10., 0.), (0., 10.), (10., 10.), (5., 5.)]
    infos = ['a', 'b', 'c', 'd', 'e']
    orig = dict(zip(pts, infos))
    print(orig)
    sorted_pts = kdsort(decorate(pts))
    print(sorted_pts)
    new2old = translation_table(sorted_pts)
    print(new2old)

    dt = Triangulation()
    incremental = KDOrderPointInserter(dt)
    incremental.insert(sorted_pts)

    for new_index, v in enumerate(dt.vertices):
        v.info = infos[new2old[new_index]]

    lst = [((v.x, v.y), v.info) for v in dt.vertices]
    print(orig)
    print(dict(lst))
    assert orig == dict(lst)


def test_empty_triangulation():
    print('empty triangulation')
    logging.basicConfig(level=logging.DEBUG)
    dt = Triangulation()
    incremental = KDOrderPointInserter(dt)
    incremental.initialize([(0, 0)])
    print("done.")
    with open("/tmp/all_tris.wkt", "w") as fh:
        print(len(dt.triangles))
        output_triangles(dt.triangles, fh)
    with open("/tmp/all_vertices.wkt", "w") as fh:
        output_vertices(dt.vertices, fh)


def test_triangulate():
    logging.basicConfig(level=logging.DEBUG)
    from tri.delaunay.helpers import ToPointsAndSegments
#     polygon = [
#         [(256.0, 760.0), (518.0, 760.0), (518.0, 630.0), (674.0, 630.0), (674.0, 239.0), (673.0, 239.0), (127.0, 239.0), (127.0, 240.0), (126.0, 240.0), (126.0, 513.0), (127.0, 513.0), (127.0, 514.0), (126.0, 514.0), (126.0, 630.0), (255.0, 630.0), (256.0, 630.0), (256.0, 760.0)],
#         [(128.0, 629.0), (128.0, 423.0), (270.0, 423.0), (270.0, 422.0), (271.0, 422.0), (271.0, 240.0), (672.0, 240.0), (672.0, 629.0), (128.0, 629.0)],
#         [(258.0, 759.0), (258.0, 631.0), (516.0, 631.0), (516.0, 759.0), (258.0, 759.0)],
#         [(128.0, 421.0), (128.0, 240.0), (269.0, 240.0), (269.0, 421.0), (128.0, 421.0)]
#     ]
    ln = [(0, 0), (2, 0), (3, 0), (4, 15), (5, 0), (7, 0), (8, 15), (9, 0), (12, 0), (12, 50), (11, 50), (10, 35), (9, 50), (7, 50), (6, 35), (5, 50), (3, 50), (2, 35), (1, 50), (0, 50), (0, 0)]

    helper = ToPointsAndSegments()
#     helper.add_polygon(polygon)
    helper.add_linestring(ln)

    triangulate(helper.points, None, helper.segments, output=True)

#     test_sort_keeps_connection_info()

#     test_sort_keeps_connection_info()
#     pts = [(0, 0), (1, 1), (1, 0), (7, 0.5)]
#     segments = [(0, 3)]
# #     from tri.delaunay.helpers import random_circle_vertices
# #     pts = set(random_circle_vertices(15000))
#     infos = ["{}-info".format(i) for (i, pt) in enumerate(pts)]
#     triangulate(pts, infos, segments, output=True)

def test_triangulate_and_remove1pt():
    pts = [(0, 0), (-15, 9)]
#    pts = [(0, 0)]

    dt = Triangulation()
    incremental = KDOrderPointInserter(dt)

    sorted_pts = kdsort(decorate(pts))

    incremental.insert(sorted_pts)

    with open("/tmp/all_tris.wkt", "w") as fh:
        print(len(dt.triangles))
        output_triangles(dt.triangles, fh)

    with open("/tmp/all_vertices.wkt", "w") as fh:
        print(len(dt.vertices))
        output_vertices(dt.vertices, fh)

    pre_len_v = len(dt.vertices)
    pre_len_t = len(dt.triangles)
    incremental.remove( (0, 0), dt.triangles[0] )
    assert len(dt.vertices) + 1 == pre_len_v
    assert len(dt.triangles) + 2 == pre_len_t

    with open("/tmp/all_tris.wkt", "w") as fh:
        print(len(dt.triangles))
        output_triangles(dt.triangles, fh)

    with open("/tmp/all_vertices.wkt", "w") as fh:
        print(len(dt.vertices))
        output_vertices(dt.vertices, fh)


def test_triangulate_and_remove1pt_circle():
    # -- construct points that are on unit circle
    from math import sin, cos, pi
    ct = 1000
    alpha = 2*pi / ct
    pts = [(0, 0)]
    for i in range(ct):
        pt = (sin(i * alpha), cos(i * alpha))
        pts.append(pt)

    # -- insert
    dt = Triangulation()
    incremental = KDOrderPointInserter(dt)
    sorted_pts = kdsort(decorate(pts))
    incremental.insert(sorted_pts)

    # -- remove
    pre_len_v = len(dt.vertices)
    pre_len_t = len(dt.triangles)
    incremental.remove( (0, 0), dt.triangles[0] )
    assert len(dt.vertices) + 1 == pre_len_v
    assert len(dt.triangles) + 2 == pre_len_t

    # -- output
    with open("/tmp/all_tris.wkt", "w") as fh:
        print(len(dt.triangles))
        output_triangles(dt.triangles, fh)
    with open("/tmp/all_vertices.wkt", "w") as fh:
        print(len(dt.vertices))
        output_vertices(dt.vertices, fh)



if __name__ == "__main__":
#    test_triangulate_and_remove1pt_circle()
#    test_empty_triangulation()
#    test_triangulate()
    test_sort_keeps_connection_info()
