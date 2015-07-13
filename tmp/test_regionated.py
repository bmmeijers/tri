from tri import ToPointsAndSegments, triangulate
from tri.delaunay import output_triangles, TriangleIterator
from tri.delaunay import RegionatedTriangleIterator
class Foo(object):
    def __init__(self):
        pass
    #@profile
    def main(self):
        pts_segs = ToPointsAndSegments()
        pts_segs.add_polygon([[(0,0), (10,0), (5,10), (0,0)],
                              [(0,0), (8,2), (6,4), (5,7), (0,0)]
                              ],
                             )
        #pts_segs.add_polygon([[(10,0), (15,10), (5,10), (10,0)],
                              #[(2,2), (8,2), (6,4), (5,7), (2,2)]
                              #],
                             #)

        dt = triangulate(pts_segs.points, pts_segs.infos, pts_segs.segments)

        with open("/tmp/alltris.wkt", "w") as fh:
            output_triangles([t for t in TriangleIterator(dt)], fh)

        with open("/tmp/path.wkt", "w") as fh:
            fh.write("i;group;depth;wkt\n")
            it = RegionatedTriangleIterator(dt)
            for i, (g, d, t) in enumerate(it, start = 1):
                fh.write("{0};{1};{2};{3}\n".format(i, g, d, t))

        with open("/tmp/later.wkt", "w") as fh:
            fh.write("wkt\n")
            for t in it.later:
                fh.write("{0}\n".format(t))

Foo().main()
