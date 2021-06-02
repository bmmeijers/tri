#from simplegeom.geometry import LineString, Polygon, LinearRing

#from topomap.topomap import TopoMap
#from topomap.loopfactory import find_loops

#from tri.delaunay import ToPointsAndSegments, triangulate
##from tri.delaunay import output_triangles,
#from tri.delaunay.iter import TriangleIterator, InteriorTriangleIterator
#from tri.delaunay.iter import RegionatedTriangleIterator, ccw, StarEdgeIterator

from tri.delaunay.tds import ccw
from tri.delaunay.iter import StarEdgeIterator

# from splitarea.flagging import MidpointHarvester
# from splitarea.densify import densify

#
# FIXME:
# - With 90 degree triangle this will not work: 2 in-centers at exactly same place
# - Happily converts a constrained Triangulation as well: can lead to problems!
# - Should it embed the construction of the TopoMap object in the transformer?
# - Does not deal with SRID for the Topomap generated yet
#

def dist_squared(orig, dest):
    dx = orig[0] - dest[0]
    dy = orig[1] - dest[1]
    return (dx * dx + dy * dy) #** 0.5


class VoronoiTransformer(object):
    """Class to transform a Delaunay triangulation into a Voronoi diagram

    The class generates a series of segments, together with information how
    these should be glued together to the Voronoi diagram 
    (start node id, end node id, left face id, right face id)
    """

    def __init__(self, triangulation):
        self.triangulation = triangulation

    def transform(self):
        """Calculate center of circumscribed circles for all triangles
        and generate a line segment from one triangle to its neighbours 
        (this happens only once for every pair).
        """
#        tri_keys = {}
#        vtx_keys = {}
#        for key, t in enumerate(self.triangulation.triangles, start = 1):
#            tri_keys[t] = key
#        for key, v in enumerate(self.triangulation.vertices, start = 1):
#            vtx_keys[v] = key
        self._transform_centers()
        self._transform_segments()

    def _transform_centers(self):
        self.centers = {}
        for t in self.triangulation.triangles:
            if t.is_finite:
                self.centers[id(t)] = self.circumcenter(t) #tri_keys[t]

    def _transform_segments(self):
        segments = []
        for t in self.triangulation.triangles:
            for side, n in enumerate(t.neighbours):
#                if n is not None and \
#                    n is not self.triangulation.external and \
#                        id(t) < id(n):
                if n.is_finite and \
                        id(t) < id(n):
                        #tri_keys[t] < tri_keys[n]:
                    start, end = id(t), id(n) #tri_keys[t], tri_keys[n]
                    # dependent on whether this is a finite or infinite vertex
                    # we set the left / right face pointer
                    left_vertex = t.vertices[ccw(ccw(side))]
                    if left_vertex.is_finite:
                        left = left_vertex #vtx_keys[left_vertex]
                    else:
                        left = None
                    right_vertex = t.vertices[ccw(side)]
                    if right_vertex.is_finite:
                        right = right_vertex #vtx_keys[right_vertex]
                    else:
                        right = None
                    segments.append((start, end, left, right))
        self.segments = segments

    def close_centers(self):
        close = {}
        for v in self.triangulation.vertices:
            # subsequent triangles that have their circumcenter on the same
            # point/are very close could mess with the topology
            star = (edge.triangle for edge in StarEdgeIterator(v))
            ring = [id(triangle) for triangle in star]
            other = ring[1:] + [ring[0]]
            pairs = zip(ring, other)
            for pair in pairs:
                a, b = pair[0], pair[1]
                if a < b:
                    d = dist_squared(self.centers[a], self.centers[b])
                    if d < 1e-10:
                        if a not in close:
                            close[a] = {b: d}
                        else:
                            close[a][b] = d
                        if b not in close:
                            close[b] = {a: d}
                        else:
                            close[b][a] = d
        return close

    def _transform_cells(self):
        cells = []
        for v in self.triangulation.vertices:
            # subsequent triangles that have their circumcenter on the same
            # point/are very close could mess with the topology
            star = (edge.triangle for edge in StarEdgeIterator(v))
            ring = [self.centers[id(triangle)] for triangle in star if triangle.is_finite]
#            poly = Polygon(LinearRing(ring + [ring[0]]))
            poly = [ring + [ring[0]]]
            cells.append(poly)
        return cells

#    def incenter(self, t):
#        """Returns for a Triangle *t* the coordinates of the in-center
#        """
#        p0, p1, p2, = t.vertices
#        # subtract p0 for more precision
#        ax, ay, bx, by, cx, cy, = 0., 0., p1.x - p0.x, p1.y - p0.y, p2.x - p0.x, p2.y - p0.y
#        a2 = pow(ax, 2) + pow(ay, 2)
#        b2 = pow(bx, 2) + pow(by, 2)
#        c2 = pow(cx, 2) + pow(cy, 2)
#        ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by))
#        uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax))
#        D = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
#        ux /= D
#        uy /= D
#        # add p0
#        ux += p0.x
#        uy += p0.y
#        return (ux, uy)


    def circumcenter(self, t):
        """Returns for a Triangle *t* the coordinates of its circumcenter
        """
        p0, p1, p2, = t.vertices
        ax, ay, bx, by, cx, cy = p0.x, p0.y, p1.x, p1.y, p2.x, p2.y
        bx -= ax
        by -= ay
        cx -= ax
        cy -= ay

        bl = bx * bx + by * by
        cl = cx * cx + cy * cy

        d = bx * cy - by * cx

        x = (cy * bl - by * cl) * 0.5 / d
        y = (bx * cl - cx * bl) * 0.5 / d

        return (ax + x, ay + y)


def main():
    import json
#     pts_segs = ToPointsAndSegments()
#     pts_segs.add_polygon([[(0,0), (10,0), (5,10), (0,0)],
#                           #[(2,2), (8,2), (6,4), (5,7), (2,2)]
#                           ],
#                          )
#     pts_segs.add_polygon([[(10,0), (15,10), (5,10), (10,0)],
#                           #[(2,2), (8,2), (6,4), (5,7), (2,2)]
#                           ],
#                          )

    # FIXME: does not work with this dataset yet, as the vertex density is not
    # high enough: should add more vertices (densify)
    with open('/home/martijn/workspace/splitarea/data/sandro/poly.geojson') as fh:
        c = json.loads(fh.read())
    conv = ToPointsAndSegments()
    poly = c['features'][0]['geometry']['coordinates']
    rings = []
    for ring in poly:
        rings.append(
#                      densify(
                     [tuple(pt) for pt in ring]
#                      , 5)
                     )
    del poly
    conv.add_polygon(rings)
    dt = triangulate(conv.points, conv.infos, conv.segments)

    trafo = VoronoiTransformer(dt)
    trafo.transform()

    with open("/tmp/vroni.wkt", "w") as fh:
        fh.write("wkt;start;end;left;right\n")
        for (start, end, lft, rgt) in trafo.segments:
            fh.write("LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]});{2};{3};{4};{5}\n".format(trafo.centers[start], trafo.centers[end], start, end, lft, rgt))

    # FIXME: this should be part of the VoronoiTransformer !
    tm = TopoMap()
    for i, (start, end, lft, rgt) in enumerate(trafo.segments, start = 1):
        tm.add_edge(i, 
                    start, end, 
                    lft, rgt, 
                    LineString([trafo.centers[start], trafo.centers[end]]))
    find_loops(tm)
    with open("/tmp/geom.wkt", "w") as fh:
        fh.write("wkt\n")
        for face in tm.faces.itervalues():
            try:
                fh.write("{0}\n".format(face.multigeometry()[0]))
            except:
                pass

#     visitor = MidpointHarvester([t for t in InteriorTriangleIterator(dt)])
#     visitor.skeleton_segments()
#     with open("/tmp/skel.wkt", "w") as fh:
#         fh.write("wkt\n")
#         for seg in visitor.segments:
#             fh.write("LINESTRING({0[0].x} {0[0].y}, {0[1].x} {0[1].y})\n".format(seg))

    with open("/tmp/inside.wkt", "w") as fh:
        output_triangles([t for t in TriangleIterator(dt)], fh)

class DisjointSet(object):
    """
    Taken from: https://stackoverflow.com/a/3067672
    """

    def __init__(self):
        self.leader = {} # maps a member to the group's leader
        self.group = {}  # maps a group leader to the group (which is a set)

    def add(self, a, b):
        leadera = self.leader.get(a)
        leaderb = self.leader.get(b)
        if leadera is not None:
            if leaderb is not None:
                if leadera == leaderb:
                    return # nothing to do
                groupa = self.group[leadera]
                groupb = self.group[leaderb]
                if len(groupa) < len(groupb):
                    # swap, add to largest group
                    a, leadera, groupa, b, leaderb, groupb = b, leaderb, groupb, a, leadera, groupa
                groupa |= groupb
                del self.group[leaderb]
                for k in groupb:
                    self.leader[k] = leadera
            else:
                self.group[leadera].add(b)
                self.leader[b] = leadera
        else:
            if leaderb is not None:
                self.group[leaderb].add(a)
                self.leader[a] = leaderb
            else:
                self.leader[a] = self.leader[b] = a
                self.group[a] = set([a, b])

def transform(dt):
    trafo = VoronoiTransformer(dt)
    trafo._transform_centers() # calculate circumcenter of all triangles
    # FIXME: this is naive: segments can be unneeded, if there are triangles that have the same circumcentre, or for which the center is exactly the same (unlikely under floating point precision -> angle calculation is also not possible then)...
    trafo._transform_segments()
    #with open("/tmp/vroni_segments.wkt", "w") as fh:
    #    fh.write("id;wkt;start;end;left;right;dist\n")
    #    for i, (start, end, lft, rgt) in enumerate(trafo.segments, start = 1):
    #        fh.write("{6};LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]});{2};{3};{4};{5};{7}\n".format(trafo.centers[start], trafo.centers[end], start, end, lft, rgt, i, dist(trafo.centers[start], trafo.centers[end])))
    # by iterating over all vertices and all the triangles in their stars we 
    # obtain the voronoi cells (polygons)
    cells = trafo._transform_cells()
    with open("/tmp/vroni_cells.wkt", "w") as fh:
        fh.write("wkt\n")
        for cell in cells:
            fh.write("{}\n".format(cell))

    close_pairs = trafo.close_centers()

    ds = DisjointSet()
    for a in close_pairs:
        for b in close_pairs[a]:
            ds.add(a, b)
    
    # print ds.leader
    # print ds.group
    
    
    # trafo.transform()
    
    
    
    #with open("/tmp/allvertices.wkt", "w") as fh:
    #    tri.delaunay.output_vertices(dt.vertices, fh)
    
    with open("/tmp/vroni_nodes.wkt", "w") as fh:
        fh.write("id;wkt\n")
        for i, geom in trafo.centers.iteritems():
            txt = "{0};POINT({1[0]} {1[1]})\n".format(i, geom)
            fh.write(txt)
    
    #def dist(orig, dest):
    #    dx = orig[0] - dest[0]
    #    dy = orig[1] - dest[1]
    #    return (dx * dx + dy * dy) ** 0.5
    
    with open("/tmp/vroni_segments_culled.wkt", "w") as fh:
        fh.write("id;wkt;start;end;left;right;dist\n")
        for i, (start, end, lft, rgt) in enumerate(trafo.segments, start = 1):
            start = ds.leader.get(start, start)
            end = ds.leader.get(end, end)
            # skip segments that do not have length
            if start != end:
                fh.write("{6};LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]});{2};{3};{4};{5};{7}\n".format(trafo.centers[start], trafo.centers[end], start, end, lft, rgt, i, dist(trafo.centers[start], trafo.centers[end])))
    
    from topomap.topomap import TopoMap
    from topomap.loopfactory import find_loops
    from simplegeom.geometry import LineString
    tm = TopoMap(universe_id=0)
    
    for v in dt.vertices:
        if v.is_finite:
            tm.add_face(id(v), attrs = v.info)
    
    
    for i, (start, end, lft, rgt) in enumerate(trafo.segments, start = 1):
        start = ds.leader.get(start, start)
        end = ds.leader.get(end, end)
        if start != end:
            tm.add_edge(i, 
                        start, end, 
                        lft, rgt, 
                        LineString([trafo.centers[start], trafo.centers[end]]))
    find_loops(tm)
    with open("/tmp/geom.wkt", "w") as fh:
        fh.write("wkt;segment;height\n")
        for face in tm.faces.itervalues():
            if not face.unbounded:
                fh.write("{0};{1[0]};{1[1]}\n".format(face.multigeometry()[0], face.attrs))


if __name__ == "__main__":
    main()
