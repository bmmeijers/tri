'''
Created on Nov 13, 2018

@author: martijn
'''
from random import shuffle
import logging
from datetime import datetime
from itertools import chain
from random import randint

from tri.delaunay.preds import orient2d, incircle
from tri.delaunay.iter import StarEdgeIterator
from tri.delaunay.tds import Triangle, Edge, ccw
from tri.delaunay.inout import output_triangles

# -----------------------------------------------------------------------------
# Constraints
#     The algorithm is described in:
#         Fast Segment Insertion and
#         Incremental Construction of Constrained Delaunay Triangulations
#         Jonathan Richard Shewchuk and Brielin C. Brown
#
#     Available from:
#         http://www.cs.berkeley.edu/~jrs/papers/inccdt.pdf
#
# @article{Shewchuk2015,
#   doi = {10.1016/j.comgeo.2015.04.006},
#   url = {https://doi.org/10.1016/j.comgeo.2015.04.006},
#   year = {2015},
#   month = sep,
#   publisher = {Elsevier {BV}},
#   volume = {48},
#   number = {8},
#   pages = {554--574},
#   author = {Jonathan Richard Shewchuk and Brielin C. Brown},
#   title = {Fast segment insertion and incremental construction of constrained Delaunay triangulations},
#   journal = {Computational Geometry}
# }


def triangle_overlaps_ray(vertex, towards):
    """Returns the triangle that overlaps the ray.
    In case there are multiple candidates,
    then the triangle with the right
    leg overlapping the ray is returned.
    It's a ValueError if no or multiple candidates are found.
    """
    candidates = []
    for edge in StarEdgeIterator(vertex):
        # if edge.isFinite:
        start, end = edge.segment
        # start: turns ccw
        # end: turns cw
        ostart = orient2d(start, towards, vertex)
        oend = orient2d(end, towards, vertex)
        if ostart >= 0 and oend <= 0:
            candidates.append((edge, ostart, oend))
    # the default, exactly one candidate
    if len(candidates) == 1:
        return candidates[0][0]
    # no candidates found,
    # this would be the case if towards lies outside
    # currently triangulated convex hull
    elif len(candidates) == 0:
        raise ValueError(
            "No overlap found (towards outside triangulated convex hull?)")
    # the ray overlaps the legs of multiple triangles
    # only return the triangle for which the right leg overlaps with the ray
    # it is an error if there is not exactly one candidate that we can return
    else:
        ostartct = 0
        candidateIdx = None
        for i, (edge, ostart, oend) in enumerate(candidates):
            if ostart == 0:
                ostartct += 1
                candidateIdx = i
        if ostartct != 1 or candidateIdx is None:
            for i, (edge, ostart, oend) in enumerate(candidates):
                print((ostart, oend))
            raise ValueError("Incorrect number of triangles found")
        return candidates[candidateIdx][0]


def mark_cavity(P, Q, triangles):
    """Returns two lists: Edges above and below the list of triangles.
    These lists are sorted clockwise around the triangles
    (this is needed for CavityCDT).
    """
    # From a list of triangles make two lists of edges:
    # above and below...
    # It is made sure that the edges that are put
    # here are forming a polyline
    # that runs *clockwise* around the cavity
    # - to output the cavity:
#    with open('/tmp/cavity.wkt', 'w') as fh:
#         output_triangles(triangles, fh)
    assert len(triangles) != 0
    above = []
    below = []
    if len(triangles) == 1:
        t = triangles[0]
        pidx = t.vertices.index(P)
        lidx = (pidx + 1) % 3
        ridx = (pidx + 2) % 3
        lv = t.vertices[lidx]
        # r = t.vertices[ridx]
#         print "p", P
#         print "q", Q
#         print "lv", lv
#         print "r", r
        assert lv is Q
#         if lv is Q:
#             print "L IS Q"
        below = []
        for i in (ridx,):
            n = t.neighbours[i]
            b = Edge(n, n.neighbours.index(t))
            # b = Edge(n, common(t, i, n))
            below.append(b)
        above = []
        for i in (lidx, pidx,):
            n = t.neighbours[i]
            # b = Edge(n, common(t, i, n))
            b = Edge(n, n.neighbours.index(t))
            above.append(b)
#             below = [Edge(t.getNeighbour(pidx), t.getOppositeSide(pidx))]
#             above = [Edge(t.getNeighbour(ridx), t.getOppositeSide(ridx)),
#                      Edge(t.getNeighbour(lidx), t.getOppositeSide(lidx))]

#         elif r is Q:
#             print "R IS Q"
#             below = []
#             for i in (pidx, lidx,):
#                 n = t.neighbours[i]
#                 b = Edge(n, common(t, i, n))
#                 below.append(b)
#             above = []
#             for i in (ridx,):
#                 n = t.neighbours[i]
#                 b = Edge(n, common(t, i, n))
#                 above.append(b)
#             above = [Edge(t.getNeighbour(ridx), t.getOppositeSide(ridx))]
#             below = [Edge(t.getNeighbour(lidx), t.getOppositeSide(lidx)),
#                      Edge(t.getNeighbour(pidx), t.getOppositeSide(pidx))
#         else:
#             raise ValueError("Combinations exhausted")
    else:
        # precondition here is that triangles their legs
        # do NOT overlap with the segment that goes
        # from P -> Q
        # thus: left and right orientation cannot both be 0
        # -> this is an error
        for t in triangles:
            for side in range(3):
                # FIXME:
                # skip neighbouring side if it is none?
#                 if t.neighbours[side] is None:
#                     continue
                # FIXME: Do we add None into above/below???
                # if that is the case, we should still know which 2 infinite
                # vertices are between here...

                edge = Edge(t, side)
                R, L = edge.segment
                left = orient2d(L, Q, P)
                right = orient2d(R, Q, P)
                # in case both are 0 ... not allowed
                if (left == 0 and right == 0):
                    raise ValueError("Overlapping triangle leg found,"
                                     " not allowed")
                n = t.neighbours[side]
                e = Edge(n, n.neighbours.index(t))
                if left >= 0 and right >= 0:
                    below.append(e)
                elif right <= 0 and left <= 0:
                    above.append(e)

        below.reverse()
    return above, below


def straight_walk(P, Q):
    """Obtain the list of triangles that overlap
    the line segment that goes from Vertex P to Q.

    Note that P and Q must be Vertex objects that are in the Triangulation
    already.

    Raises a ValueError when either a Constrained edge is crossed in the
    interior of the line segment or when another Vertex lies on the
    segment.
    """
    edge = triangle_overlaps_ray(P, Q)
    t = edge.triangle
    side = edge.side
    R, L = edge.segment
    out = [t]
    if Q in t.vertices:
        # we do not need to go into walking mode if we found
        # the exact triangle with the end point already
        return out
    # perform walk
    # pos = t.vertices.index(R)

    # from end via right to left makes right turn (negative)
    # if line is collinear with end point then orientation becomes 0

    # FIXME:
    # The way that we now stop the rotation around the vertex
    # does that make a problem here --> we can get either the lower
    # or the upper triangle, this depends on the arbitrary start triangle
    while orient2d(Q, R, L) < 0.:
        # check if we do not prematurely have a orientation of 0
        # at either side, which means that we collide a vertex
        if (L is not Q and orient2d(L, P, Q) == 0) or \
                (R is not Q and orient2d(R, P, Q) == 0):
            raise ValueError("Unwanted vertex collision detected - inserting: {} -> {} | crossing: {} -> {}". format(P, Q, R, L))

        # based on the position of R take next triangle
        # FIXME:
        # TEST THIS: check here if taking the neighbour does not take
        # place over a constrained side of the triangle -->
        # raise ValueError("Unwanted constrained segment collision detected")
#         if triangle.getEdgeType(side):
#             raise TopologyViolationError("Unwanted"
#                        " constrained segment collision detected")
        if t.constrained[side]:
            raise ValueError("Unwanted constrained segment collision detected - inserting: {} -> {} | crossing: {} -> {}". format(P, Q, R, L))
        t = t.neighbours[side]
        out.append(t)

        side = t.vertices.index(R)
        S = t.vertices[ccw(side)]
        ori = orient2d(S, Q, P)
        #
        if ori < 0:
            L = S
            side = ccw(side+1)
        else:
            R = S
        # check if we do not prematurely have a orientation of 0
        # at either side, which means that we collide a vertex
        if (L is not Q and orient2d(L, P, Q) == 0) or \
                (R is not Q and orient2d(R, P, Q) == 0):
            raise ValueError("Unwanted vertex collision detected - inserting: {} -> {} | crossing: {} -> {}". format(P, Q, R, L))

    return out


def permute(a, b, c):
    """Permutation of the triangle vertex indices from lowest to highest,
    i.e. a < b < c

    This order makes sure that a triangle is always addressed in the same way

    Used in CavityCDT.
    """
    return tuple(sorted([a, b, c]))


class ConstraintInserter(object):
    """Constraint Inserter

    Insert segments into a Delaunay Triangulation.
    """

    def __init__(self, triangulation):
        self.dt = triangulation

    def insert(self, segments):
        """Insert constraints into dt

        Parameter: segments - list of 2-tuples, with coordinate indexes
        """
        for j, segment in enumerate(segments):
            p = self.dt.vertices[segment[0]]
            q = self.dt.vertices[segment[1]]
            try:
                self.insert_constraint(p, q)
            except Exception:
                logging.debug("skipping {} {}".format(p, q))
                raise

            if (j % 10000) == 0:
                logging.debug(" " + str(datetime.now()) + str(j))
        self.remove_empty_triangles()

    def remove_empty_triangles(self):
        """Removes empty triangles (not pointing to any vertex) by filtering
        the triangles that have one of its vertex members set
        """
        new = [x for x in self.dt.triangles if not(x.vertices[0] is None or x.vertices[1]
                                   is None or x.vertices[2] is None)]
        self.dt.triangles = new

    def insert_constraint(self, P, Q):
        """Insert constraint into dt.

        It leaves the triangles that are removed inside the cavity of
        the constraint inserted in the triangles array
        """
        logging.debug(" constraint LINESTRING({}, {})".format(P, Q))
        if P is Q:
            logging.warning("Equal points found while "
                            "inserting constraint: {} {} -- skipped insertion".format(P, Q))
            return
#             raise ValueError(
#                 "Equal points found "
#                 "inserting constraint: {} {}".format(P, Q))
        cavity = straight_walk(P, Q)
        above, below = mark_cavity(P, Q, cavity)
        # change triangle pointers of vertices around the cavity to point to
        # triangles that lie outside the cavity (as these will be removed
        # later)
        for edge in chain(above, below):
            for vertex in edge.segment:
                vertex.triangle = edge.triangle
        # Re-triangulate upper half
        cavA = CavityCDT(self.dt, above)
        A = cavA.edge
        # Re-triangulate bottom half
        cavB = CavityCDT(self.dt, below)
        B = cavB.edge
        # link up the two triangles at both sides of the segment
        A.triangle.neighbours[A.side] = B.triangle
        B.triangle.neighbours[B.side] = A.triangle
        # constrained edges
        A.triangle.constrained[A.side] = True
        B.triangle.constrained[B.side] = True
        for t in cavity:
            t.vertices = [None, None, None]
            t.neighbours = [None, None, None]


class CavityCDT(object):
    """Class to triangulate an `evacuated' cavity adjacent to a constraint
    """

    def __init__(self,
                 triangulation,
                 cavity_edges):
        """
        dt - the triangulation data structure
        cavity_edges - the edges that bound the cavity
        in *CLOCKWISE* order
        around the cavity. Note: these edges do not include the segment
        to be inserted.
        """
        # WARNING: The ordering of vertices
        # around the cavity is important to function correctly!
        self.vertices = []
        self.edge = None
        self.dt = triangulation

        # If we found exactly one cavity edge, there is no
        # area between ray and cavity polygon. Hence this edge
        # should be the one that will be linked to (after that we've
        # set the type of this edge to constrained).
        if len(cavity_edges) == 1:
            edge = cavity_edges[0]
            # will be carried out by caller
            # edge.triangle.setEdgeType(edge.side, True)
            self.edge = edge
            return
        self._preprocess(cavity_edges)
        self._retriangulate()
        self._push_back_triangles()

    def _preprocess(self, cavity_edges):
        """Set up data structures needed for the re-triangulate part of the
        algorithm.
        """
        self.constraints = set()
        for i, edge in enumerate(cavity_edges):
            xx, yy = edge.segment
            # Both directions are needed, as this is used
            # for internal dangling edges inside the cavity,
            # which are traversed both sides.
            self.constraints.add((id(xx), id(yy)))
            self.constraints.add((id(yy), id(xx)))
            if i:
                self.vertices.append(yy)
            else:
                self.vertices.extend([xx, yy])
        # Make the vertices list COUNTERCLOCKWISE here
        # The algorithm depends on this orientation!
        self.vertices.reverse()
        self.surroundings = {}
        for i, edge in enumerate(cavity_edges):
            s = edge.segment
            self.surroundings[id(s[0]), id(s[1])] = edge
        # Make a "linked list" of polygon vertices
        self.next = {}
        self.prev = {}
        # Relative size of distances to the segment
        self.distance = {}
        # Adjacency: third point of a triangle by given oriented side
        self.adjacency = {}
        # Set of resulting triangles (vertex indices)
        self.triangles = set()
        # Initialization for the algorithm
        m = len(self.vertices)
        # Make random permutation of point indices
        self.pi = list(range(1, m - 1))
        # Randomize processing order
        shuffle(self.pi)
        # Link all vertices in a circular list that
        # describes the polygon outline of the cavity
        for i in range(m):
            self.next[i] = (i + 1) % m
            self.prev[i] = (i - 1) % m
            # Distance to the segment from [0-m]
            self.distance[i] = orient2d(self.vertices[0],
                                        self.vertices[i],
                                        self.vertices[m-1])

    def _retriangulate(self):
        """Re-triangulate the cavity, the result is a collection of
        triangles that can be pushed back into the original DT data
        structure that replaces the old triangles inside the cavity.
        """
        # Now determine how to `remove' vertices
        # from the outline in random order
        #
        # Go over pi from back to start; quit at *second* item in pi
        # This determines order of removal of vertices from
        # the cavity outline polygon
        m = len(self.vertices)
        for i in range(len(self.pi) - 1, 0, -1):
            while self.distance[self.pi[i]] < \
                        self.distance[self.prev[self.pi[i]]] and \
                        self.distance[self.pi[i]] < \
                        self.distance[self.next[self.pi[i]]]:
                # FIXME: is j correct ??? should i be i + 1 ?
                j = randint(0, i)
                self.pi[i], self.pi[j] = self.pi[j], self.pi[i]
            # take a vertex out of the circular list
            self.next[self.prev[self.pi[i]]] = self.next[self.pi[i]]
            self.prev[self.next[self.pi[i]]] = self.prev[self.pi[i]]
        # add an initial triangle
        self._add_triangle(0, self.pi[0], m-1)
        # Work through the settled order of vertex additions
        # Now in forward direction, keep adding points until all points
        # are added to the dt of this part of the cavity
        for i in range(1, len(self.pi)):
            a = self.pi[i]
            b, c = self.next[a], self.prev[a]
            self._insert_vertex(a, b, c)

    def _push_back_triangles(self):
        """Make new triangles that are inserted in the data structure
        and that are linked up properly with each other and the surroundings.
        """

        # First make new triangles
        newtris = {}
        for three in self.triangles:
            a, b, c, = three
            # Index triangle by temporary sorted list of vertex indices
            #
            # By indexing triangles this way, we
            T = Triangle(self.vertices[a],
                         self.vertices[b],
                         self.vertices[c])
            newtris[permute(id(T.vertices[0]), id(
                T.vertices[1]), id(T.vertices[2]))] = T
        for x in newtris.values():
            assert orient2d(x.vertices[0], x.vertices[1], x.vertices[2]) > 0
        # Translate adjacency table to new indices of triangles
        # Note that vertices that are used twice (because of dangling edge
        # in the cavity) will get the same identifier again
        # (while previously they would have different id's).
        adj = {}
        for (f, t), v in self.adjacency.items():
            adj[id(self.vertices[f]), id(self.vertices[t])] = id(
                self.vertices[v])
        # Link all the 3 sides of the new triangles properly
        for T in newtris.values():
            for i in range(3):
                segment = Edge(T, i).segment
                side = (id(segment[1]), id(segment[0]))
                constrained = False
                # The side is adjacent to another new triangle
                # The neighbouring triangle at this side will be linked later
                # In case this is a dangling segment we constrain the segment
                if side in adj:
                    neighbour = newtris[permute(side[0], side[1], adj[side])]
                    if side in self.constraints:
                        constrained = True
                # the side is adjacent to an exterior triangle
                # that lies outside the cavity and will
                # remain after the re-dt
                # therefore also change the neighbour of this triangle
                elif side in self.surroundings:
                    neighbour_side = self.surroundings[side].side
                    neighbour = self.surroundings[side].triangle
                    neighbour.neighbours[neighbour_side] = T  # MM fixed
                    # getEdgeType(neighbour_side)
                    constrained = neighbour.constrained[neighbour_side]
                # the triangle is the bottom of the evacuated cavity
                # hence it should be linked later to the other
                # re-dt of the cavity
                else:
                    if self.edge:
                        print((self.edge.triangle, self.edge.triangle.all_infinite, self.edge.triangle.__str__()))
                    assert self.edge is None
                    neighbour = None
                    self.edge = Edge(T, i)
                T.neighbours[i] = neighbour  # setNeighbour(i, neighbour)
                # T.setEdgeType(i, constrained)
                T.constrained[i] = constrained
            # Append the new triangles to the triangle list of the
            # dt
            self.dt.triangles.append(T)
        assert self.edge is not None

    def _insert_vertex(self, u, v, w):
        """Insert a vertex to the triangulated area,
        while keeping the area of the current polygon triangulated
        """
        x = -1
        # Find third vertex in the triangle that has edge (w, v)
        if (w, v) in self.adjacency:
            x = self.adjacency[(w, v)]
        # See if we have to remove some triangle(s) already there,
        # or that we can add just a new one
        if x != -1 and \
            (orient2d(self.vertices[u],
                      self.vertices[v],
                      self.vertices[w]) <= 0 or
             incircle(self.vertices[u],
                      self.vertices[v],
                      self.vertices[w],
                      self.vertices[x]) > 0):
            # Remove triangle (w,v,x), also from adjacency dict
            self.triangles.remove(permute(w, v, x))
            del self.adjacency[(w, v)]
            del self.adjacency[(v, x)]
            del self.adjacency[(x, w)]
            # Recurse
            self._insert_vertex(u, v, x)
            self._insert_vertex(u, x, w)
        else:
            # Add a triangle (this triangle could be removed later)
            self._add_triangle(u, v, w)

    def _add_triangle(self, a, b, c):
        """Add a triangle to the temporary set of triangles

        It is not said that a triangle that is added,
        survives until the end of the algorithm
        """
        t = permute(a, b, c)
        P = {}
        P[(a, b)] = c
        P[(b, c)] = a
        P[(c, a)] = b
        # .update() overwrites existing keys
        # (but these should not exist anyway)
        self.adjacency.update(P)
        # A triangle is stored with vertices in ordered indices
        self.triangles.add(t)
