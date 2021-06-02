'''
Created on Nov 13, 2018

@author: martijn
'''

# ------------------------------------------------------------------------------
# Iterators
#
from tri.delaunay.tds import Edge, ccw


class FiniteEdgeIterator(object):
    def __init__(self, triangulation, constraints_only=False):
        self.triangulation = triangulation
        self.constraints_only = constraints_only
        self.current_idx = 0  # this is index in the list
        self.pos = -1  # this is index in the triangle (side)

    def __iter__(self):
        return self

    def next(self):
        ret = None
        while self.current_idx < len(self.triangulation.triangles):
            triangle = self.triangulation.triangles[self.current_idx]
            # skip this triangle if it is an infinite triangle
            if not triangle.is_finite:
                self.pos = -1
                self.current_idx += 1
                continue
            self.pos += 1
            neighbour = triangle.neighbours[self.pos]
            # output edges only once:
            # inside the triangulation only the triangle
            # with lowest id its edge is output
            # along the convex hull we always output the edge
            if id(triangle) < id(neighbour) or not neighbour.is_finite:
                if self.constraints_only and triangle.constrained[self.pos]:
                    ret = Edge(triangle, self.pos)
                elif not self.constraints_only:
                    ret = Edge(triangle, self.pos)
            if self.pos == 2:
                self.pos = -1
                self.current_idx += 1
            if ret is not None:
                return ret
        else:
            raise StopIteration()

    def __next__(self):
        return self.next()


class TriangleIterator(object):
    """Iterator over all triangles that are in the triangle data structure.
    The finite_only parameter determines whether only the triangles in the
    convex hull of the point set are iterated over, or whether also infinite
    triangles are considered.

    """

    def __init__(self, triangulation, finite_only=False):
        self.triangulation = triangulation
        self.finite_only = finite_only
        self.visited = set()
        self.to_visit_stack = [self.triangulation.triangles[0]]

    def __iter__(self):
        return self

    def next(self):
        ret = None
        while self.to_visit_stack:
            triangle = self.to_visit_stack.pop()
            # determine whether we should 'emit' the triangle
            if self.finite_only and id(triangle) not in self.visited \
                    and triangle.is_finite:
                ret = triangle
            elif not self.finite_only and id(triangle) not in self.visited:
                ret = triangle
            self.visited.add(id(triangle))
            # NOTE: from an external triangle we can get
            # to a triangle in the triangulation multiple times
            for i in range(3):
                neighbour = triangle.neighbours[i]
                if neighbour is None:
                    continue
                elif id(neighbour) not in self.visited:
                    self.to_visit_stack.append(neighbour)
            if ret is not None:
                return ret
        else:
            raise StopIteration()

    def __next__(self):
        return self.next()

class ConvexHullTriangleIterator(TriangleIterator):
    """Iterator over all triangles that are in the convex hull of the
    point set (excludes infinite triangles).

    """

    def __init__(self, triangulation):
        # Actually, we are an alias for TriangleIterator
        # with finite_only set to True
        super(ConvexHullTriangleIterator, self).__init__(triangulation, True)


class InteriorTriangleIterator(object):
    """Iterator over all triangles that are enclosed by constraints

    Assumes that 1 polygon (the polygon consists of *exactly one* connected
    component) has been triangulated and the boundary is closed properly!

    Hence, pre-condition to function correctly:
    At least 1 interior triangle should exist!
    """

    def __init__(self, triangulation):
        constrained = False
        self.triangulation = triangulation
        self.visited = set()
        # start on a *INfinite* triangle
        # FIXME:
        # Potentially this is an expensive scan
        # should we carry it out once, after the triangulation algorithm is
        # finished and store the start as external
        # triangle to the triangulation object???
        # OR SHOULD WE KEEP AN EXTERNAL TRIANGLE THROUGHOUT, AS BEFORE?

        # FIXME:
        # apparently triangles[0] is not always an infinite triangle
#        start = self.triangulation.triangles[0]

        for start in self.triangulation.triangles:
             if not start.is_finite:
                 break
        assert not start.is_finite  # should be infinite triangle

        self.to_visit_stack = [start]
        # walk to an interior triangle
        while not constrained and self.to_visit_stack:
            triangle = self.to_visit_stack.pop()
            assert triangle is not None
            self.visited.add(triangle)
            # NOTE: from an external triangle we can get
            # to a triangle in the triangulation multiple times
            for i in range(3):
                constrained = triangle.constrained[i]
                neighbour = triangle.neighbours[i]
                if constrained:
                    # empty to_visit_stack
                    # start visited empty
                    self.to_visit_stack = [neighbour]
                    self.visited = set()
                    break
                if neighbour is not None and neighbour not in self.visited:
                    self.to_visit_stack.append(neighbour)

    def __iter__(self):
        return self

    def next(self):
        ret = None
        constrained = False
        while self.to_visit_stack:
            triangle = self.to_visit_stack.pop()
            if triangle not in self.visited:
                ret = triangle
            self.visited.add(triangle)
            # NOTE: from an external triangle we can get
            # to a triangle in the triangulation multiple times
            for i in range(3):
                constrained = triangle.constrained[i]
                if constrained:
                    continue
                neighbour = triangle.neighbours[i]
                if neighbour is not None and neighbour not in self.visited:
                    self.to_visit_stack.append(neighbour)
            if ret is not None:
                return ret
        else:
            raise StopIteration()

    def __next__(self):
        return self.next()


class RegionatedTriangleIterator(object):
    """Iterator over all triangles that are fenced off by constraints.
    The constraints fencing off triangles determine the regions.
    The iterator yields a tuple: (region number, depth, triangle).

    Note:

    - The region number can increase in unexpected ways, e.g. 0, 1, 476, 1440,
    ..., etc.
    - The depth gives the nesting of the regions.

    The first group is always the infinite part (at depth 0) of the domain
    around the feature (the parts of the convex hull not belonging to any
    interior part).
    """

    def __init__(self, triangulation):
        # start at the exterior
        self.triangulation = triangulation
        # set([id(self.triangulation.triangles[0])])
        # FIXME: DOES THIS HAVE A BUG? Sometimes id() is used, sometimes not???
        self.visited = set([None])
        for start in self.triangulation.triangles:
            if not start.is_finite:
                break
        else:
            raise ValueError('no infinite triangle found to start walk')
        self.to_visit_stack = [(start, 0)]
        self.later = []
        self.group = 0

    def __iter__(self):
        return self

    def next(self):
        while self.to_visit_stack or self.later:
            # visit all triangles in the exterior, subsequently visit
            # all triangles that are enclosed by a set of segments
            while self.to_visit_stack:
                triangle, depth = self.to_visit_stack.pop()
                assert triangle is not None
                if triangle in self.visited:
                    continue
                self.visited.add(triangle)
                for i in range(3):
                    constrained = triangle.constrained[i]
                    neighbour = triangle.neighbours[i]
                    if constrained and neighbour not in self.visited:
                        self.later.append((neighbour, depth + 1))
                    elif neighbour is not None and \
                            neighbour not in self.visited:
                        self.to_visit_stack.append((neighbour, depth))
                return (self.group, depth, triangle)
            # flip the next level with this
            while self.later:
                self.group += 1
                t, d = self.later.pop()
                if t not in self.visited:
                    self.to_visit_stack = [(t, d)]
                    break
        else:
            raise StopIteration()

    def __next__(self):
        return self.next()


class StarEdgeIterator(object):
    """Returns iterator over edges in the star of the vertex

    The edges are returned in counterclockwise order around the vertex.
    The triangles that the edges are associated with share the vertex
    that this iterator is constructed with.
    """

    def __init__(self, vertex):  # , finite_only = True):
        self.vertex = vertex
        self.start = vertex.triangle
        self.triangle = self.start
        self.side = ccw(self.start.vertices.index(self.vertex))
        self.done = False

    def __iter__(self):
        return self

    def next(self):
        if not self.done:
            self.triangle = self.triangle.neighbours[self.side]
            assert self.triangle is not None
            # self.visited.append(self.triangle)
            # try:
            side = self.triangle.vertices.index(self.vertex)
            # except ValueError, err:
            #    print err
            #    print [id(t) for t in self.visited]
            #    raise
            # side = (self.side + 1) % 3
            assert self.triangle.vertices[side] is self.vertex
            e = Edge(self.triangle, side)
            self.side = ccw(side)
            if self.triangle is self.start:
                self.done = True
            return e
        else:  # we are at start again
            raise StopIteration()

    def __next__(self):
        return self.next()


# class StarTriangleIterator(object):
#    """Returns iterator over triangles in the star of the vertex

#    The triangles are returned in counterclockwise order around the vertex.
#    The triangles share the vertex that this iterator is constructed with.
#    """
#    def __init__(self, vertex):
#        self.vertex = vertex
#        self.start = vertex.triangle
#        self.triangle = self.start
#        self.side = ccw(self.start.vertices.index(self.vertex))
#        self.done = False

#    def __iter__(self):
#        return self

#    def next(self):
#        if not self.done:
#            self.triangle = self.triangle.neighbours[self.side]
#            assert self.triangle is not None
#            side = self.triangle.vertices.index(self.vertex)
#            assert self.triangle.vertices[side] is self.vertex
#            self.side = ccw(side)
#            if self.triangle is self.start:
#                self.done = True
#            return self.triangle
#        else: # we are at start again
#            raise StopIteration()
