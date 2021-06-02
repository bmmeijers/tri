'''
Created on Nov 13, 2018

@author: martijn
'''
from math import sqrt, pi, cos, sin
from random import randint, random
import logging
# ------------------------------------------------------------------------------
# Generate randomized point sets (for testing purposes)
#


def random_sorted_vertices(n=10):
    """Returns a list with n random vertices
    """
    W = float(n)
    vertices = []
    for _ in range(n):
        x = randint(0, W)
        y = randint(0, W)
        x /= W
        y /= W
        vertices.append((x, y))
    vertices = list(set(vertices))
    vertices.sort()
    return vertices


def random_circle_vertices(n=10, cx=0, cy=0):
    """Returns a list with n random vertices in a circle

    Method according to:

    http://www.anderswallin.net/2009/05/uniform-random-points-in-a-circle-using-polar-coordinates/
    """
#     import fractions
#     from gmpy2 import mpfr
#     import gmpy2
#     gmpy2.get_context().precision = 53 * 4

    vertices = []
    for _ in range(n):
        r = sqrt(random())
        t = 2 * pi * random()
        x = r * cos(t)
        y = r * sin(t)
        vertices.append((x+cx, y+cy))
    vertices = list(set(vertices))
    vertices.sort()
    return vertices


class ToPointsAndSegments(object):
    """Helper class to convert a set of polygons to points and segments.
    De-dups duplicate points.
    Handles duplicate 'info' for vertices by putting the info field as a list
    to every vertex.
    """
    # FIXME: API is different because of putting info fields in list!
    def __init__(self):
        self.points = []
        self.segments = []
        self.infos = [] 
        self._points_idx = {}
        self._segments_idx = {}

    def add_polygon(self, polygon, info = None):
        """Add a polygon its points and segments to the global collection

        A polygon is a list of lists (rings), where every ring contains vertex
        objects (e.g. tuples with 2 elements).
        Important: The first and last point of a ring have to be the same 
        vertex.

        For all points in the polygon the same info will be added.
        """
        for ring in polygon:
            assert ring[0] == ring[-1]
            for pt in ring[:-1]: # skip last point of ring; should be duplicate of first
                self.add_point(pt, info)
            for start, end in zip(ring[:-1], ring[1:]):
                self.add_segment(start, end)

    def add_point(self, point, info = None):
        """Add a point and its info.
        """
        point = tuple(map(float, point))
        # -- point is not present
        if point not in self._points_idx:
            idx = len(self.points)
            self._points_idx[point] = idx
            self.points.append(point)
            # -- info given
            #    make new entry in infos list
            if info is not None:
                #self.infos.append((idx, [info]))
                self.infos.append(info)
            else:
                # no info given for this point, make empty list
                #self.infos.append((idx, []))
                #self.infos.append([])
                pass
        else:
            idx = self._points_idx[point]
            # -- add info of this point to the info list
            #if info is not None:
            #    self.infos[idx].append(info)
        return idx

    def add_linestring(self, ln, info = None):
        """Add a linestring its points.
        If info is given it is added to all the points of the line.
        """
        for pt in ln:
            self.add_point(pt, info)
        for start, end in zip(ln, ln[1:]):
            self.add_segment(start, end)

    def add_segment(self, start, end):
        """Add a segment.
        Note that points should have been added before.
        """
        start = tuple(map(float, start))
        end = tuple(map(float, end))
        start_idx, end_idx = self._points_idx[start], self._points_idx[end]
        swapped = False
        if start_idx < end_idx:
            seg = (start_idx, end_idx)
        elif start_idx > end_idx:
            seg = (end_idx, start_idx)
            swapped = True
        else:
            raise ValueError('same start as end point')
        if seg not in self._segments_idx:
            idx = len(self.segments)
            self._segments_idx[seg] = idx
            self.segments.append(seg)
        else:
            idx = self._segments_idx[seg]
        # returns signed segment index
        # (if reverse return 2-complement using ~)
        if swapped:
            return ~idx
        else:
            return idx


# class ToPointsAndSegments(object):
#     """Helper class to convert a set of polygons to points and segments.
#     De-dups duplicate points.
#     """

#     def __init__(self):
#         self.points = []
#         self.segments = []
#         self.infos = []
#         self._points_idx = {}

#     def add_polygon(self, polygon):
#         for ring in polygon:
#             logging.debug(ring)
#             self.add_linestring(ring)

#     def add_linestring(self, ln):
#         """Add a linestring its points and segments to the global collection
#         """
#         for pt in ln:
#             self.add_point(pt)
#         for start, end in zip(ln, ln[1:]):
#             self.add_segment(start, end)

#     def add_point(self, point, info=None):
#         """Add a point and its info.

#         Note that if a point already is present,
#         it is not appended nor is its info added to the infos list.
#         """
#         if point not in self._points_idx:
#             idx = len(self.points)
#             self._points_idx[point] = idx
#             self.points.append(point)
#             if info is not None:
#                 self.infos.append((idx, info))
#         else:
#             idx = self._points_idx[point]
#         return idx

#     def add_segment(self, start, end):
#         """Add a segment. Note that points should have been added before
#         """
#         self.segments.append((self._points_idx[start], self._points_idx[end]))
