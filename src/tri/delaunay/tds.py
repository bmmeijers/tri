'''
Created on Nov 13, 2018

@author: martijn
'''
from math import hypot

from tri.delaunay.preds import orient2d
# ------------------------------------------------------------------------------
# Helpers
#


# -- helper functions, could be inlined in Cythonized version
def box(points):
    """Obtain a tight fitting axis-aligned box around point set"""
    xmin = min(points, key=lambda x: x[0])[0]
    ymin = min(points, key=lambda x: x[1])[1]
    xmax = max(points, key=lambda x: x[0])[0]
    ymax = max(points, key=lambda x: x[1])[1]
    return (xmin, ymin), (xmax, ymax)


def ccw(i):
    """Get index (0, 1 or 2) increased with one (ccw)"""
    return (i + 1) % 3


def cw(i):
    """Get index (0, 1 or 2) decreased with one (cw)"""
    return (i - 1) % 3


def apex(side):
    """Given a side, give the apex of the triangle """
    return side % 3


def orig(side):
    """Given a side, give the origin of the triangle """
    return (side + 1) % 3  # ccw(side)


def dest(side):
    """Given a side, give the destination of the triangle """
    return (side - 1) % 3  # cw(side)


class Vertex(object):
    """A vertex in the triangulation.
    Can carry extra information via its info property.
    """
    __slots__ = ('x', 'y', 'info', 'triangle')

    def __init__(self, x, y, info=None):
        self.x = x
        self.y = y
        self.info = info
        self.triangle = None

    def __str__(self):
        return "{0} {1}".format(self.x, self.y)

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        else:
            raise IndexError("No such ordinate: {}".format(i))

    def __eq__(self, other):
        if other is None:
            return False
        return True if self.x == other.x and self.y == other.y else False

    def __hash__(self):
        return hash((self.x, self.y))


    def distance(self, other):
        """Cartesian distance to other point """
        # only used in triangle.__str__
        return hypot(self.x - other.x, self.y - other.y)

    def distance2(self, other):
        """Cartesian distance *squared* to other point """
        # Used for distances in random triangle close to point
        return pow(self.x - other.x, 2) + pow(self.y - other.y, 2)

    @property
    def is_finite(self):
        return True


class InfiniteVertex(Vertex):
    __slots__ = ('x', 'y', 'info', 'triangle')

    def __init__(self, x, y, info=None):
        super(InfiniteVertex, self)
        self.x = float(x)
        self.y = float(y)
        self.info = info
        self.triangle = None

    @property
    def is_finite(self):
        return False

#     @property
#     def x(self):
#         raise ValueError("Infinite vertex has no geometric embedding")
#
#     @property
#     def y(self):
#         raise ValueError("Infinite vertex has no geometric embedding")

#     def __str__(self):
#         return u"Inf Inf"


class Triangle(object):
    """Triangle for which its vertices should be oriented CCW
    """

    __slots__ = ('vertices', 'neighbours', 'constrained', 'info')

    def __init__(self, a, b, c):
        self.vertices = [a, b, c]  # orig, dest, apex -- ccw
        self.neighbours = [None] * 3
        self.constrained = [False] * 3  # FIXME: as bitmask on an integer
        self.info = None

    def __str__(self):
        """Conversion to WKT string.
        Defines a geometric embedding of the Infinite vertex
        so that the vertex lies perpendicular halfway convex hull segment
        """
        vertices = []
        for idx in range(3):
            v = self.vertices[idx]
            if v is not None:
                vertices.append(str(v))
            else:
                orig_idx, dest_idx = (idx - 1) % 3, (idx + 1) % 3
                orig, dest = self.vertices[orig_idx], self.vertices[dest_idx]
                halfway = (orig.x + dest.x) * .5, (orig.y + dest.y) * .5
#                 print(halfway)
                d = orig.distance(dest)
                dx = dest.x - orig.x
#                 print(d)
#                 print(dx)
                dx /= d
                dy = dest.y - orig.y
#                 print(dy)
                dy /= d
                dx *= d
                dy *= d
                pt_halfway = halfway[0] + dy, halfway[1] - dx
#                 print("outside", orig_idx, dest_idx, pt_halfway)
                vertices.append("{0[0]} {0[1]}".format(pt_halfway))
        vertices.append(vertices[0])
        return "POLYGON(({0}))".format(", ".join(vertices))

    @property
    def is_finite(self):
        return not any([isinstance(v, InfiniteVertex) for v in self.vertices])

    @property
    def all_infinite(self):
        return all([isinstance(v, InfiniteVertex) for v in self.vertices])

    @property
    def external(self):
        return self.vertices[2] is None

    @property
    def is_ccw(self):
        return orient2d(self.vertices[0],
                        self.vertices[1],
                        self.vertices[2]) > 0.


class Edge(object):
    """An edge is a Triangle and an integer [0, 1, 2] that indicates the
    side of the triangle to use as the Edge"""

    def __init__(self, triangle, side):
        self.triangle = triangle
        self.side = side

    def __eq__(self, other):
        return True if self.triangle is other.triangle and \
            self.side == other.side else False

    @property
    def segment(self):
        return (self.triangle.vertices[ccw(self.side)],
                self.triangle.vertices[cw(self.side)])

    @property
    def constrained(self):
        return self.triangle.constrained[self.side]

#     @property
#     def is_finite(self):
#         # FIXME:
#         # not use triangle here, but check if vertices are finite or infinite
#         return self.triangle.is_finite


class Triangulation(object):
    """Triangulation data structure"""
    # This represents the mesh

    def __init__(self):
        self.vertices = []
        self.triangles = []
#         self.external = None
        # infinite, external triangle (outside convex hull)
