'''
Created on Nov 13, 2018

@author: martijn
'''



# -----------------------------------------------------------------------------
# Test methods
#

def test_circle():
    """Test points in some clusters.
    """
    n = 15000
    vertices = random_circle_vertices(n, 0, 0)
#     vertices.extend(random_circle_vertices(n, 3, 4.5))
#     vertices.extend(random_circle_vertices(n, 4.6, 0.2))
#     vertices.extend(random_circle_vertices(n, 7, 2.5))
#     vertices.extend(random_circle_vertices(n, 5, -5))
#     vertices.extend(random_circle_vertices(n, 10, 5))
#     vertices.extend(random_circle_vertices(n, 9, -1))
#     vertices.extend(random_circle_vertices(n, 15, -5))
    dt = triangulate(vertices)
    for v in dt.vertices:
        print(v.x)

def test_incremental():
    L = random_sorted_vertices(n = 125000)
    tds = triangulate(L)
    with open("/tmp/alltris.wkt", "w") as fh:
                output_triangles([t for t in TriangleIterator(tds, 
                                                              finite_only=False)], 
                                 fh)
    with open("/tmp/allvertices.wkt", "w") as fh:
        output_vertices(tds.vertices, fh)

def test_cpo():
    # i = 1, j = 100 -> worst case, all end up as 1 point in slot
    # --> would be better to switch then to other order
    points = []
    idx = 0
    for i in range(400):
        for j in range(400):
            points.append((i, j, None, idx))
            idx += 1
    points_hcpo = hcpo(points)
    assert len(points) == len(points_hcpo)
    #print points
    # build a translation table for indices in the points list
#     index_translation = dict([(newpos, pos) for (newpos, (_, _, _, pos)) in enumerate(points_hcpo)])
    #print index_translation
#     with open("/tmp/points.txt", "w") as fh:
#         print >> fh, "i;wkt"
#         for i, pt in enumerate(points_hcpo):
#             print >> fh, i, ";POINT({0[0]} {0[1]})".format(pt)

def test_square():
    triangulate([(0.,0.), (10.,0.), (10., 10.), (0.,10.)], 
                [(0,1), (1,2), (2,3), (3,0)])



def test_poly():
    from connection import connection
    db = connection(True)

    def polygon_input(lines):
        points = []
        segments = []
        points_idx = {}
        for line in lines:
            for pt in line:
                if pt not in points_idx:
                    points_idx[pt] = len(points)
                    points.append(pt)
            for start, end in zip(line[:-1], line[1:]):
                segments.append((points_idx[start], points_idx[end]))
        return points, segments

    lines = []
    sql = 'select geometry from clc_edge where left_face_id in (45347) or right_face_id in (45347)'
    #sql = 'select geometry from clc_edge where left_face_id in (28875) or right_face_id in (28875)'
    # 45270
    sql = 'select geometry from clc_edge where left_face_id in (45270) or right_face_id in (45270)'
    for geom, in db.recordset(sql):
        lines.append(geom)
    points, segments = polygon_input(lines)
    dt = triangulate(points, segments)
    #
    if False:
        trafo = VoronoiTransformer(dt)
        trafo.transform()
        with open("/tmp/centers.wkt", "w") as fh:
            fh.write("wkt\n")
            for incenter in trafo.centers.itervalues():
                fh.write("POINT({0[0]} {0[1]})\n".format(incenter))
        with open("/tmp/segments.wkt", "w") as fh:
            fh.write("wkt\n")
            for segment in trafo.segments:
                # FIXME: why are some not coming through?
                try:
                    fh.write("LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})\n".format(trafo.centers[segment[0]],
                                                                             trafo.centers[segment[1]]))
                except:
                    pass
    if True:
        with open("/tmp/alltris.wkt", "w") as fh:
                    output_triangles([t for t in TriangleIterator(dt, 
                                                                  finite_only=False)], 
                                     fh)
        with open("/tmp/allvertices.wkt", "w") as fh:
            output_vertices(dt.vertices, fh)
        with open("/tmp/interiortris.wkt", "w") as fh:
                    output_triangles([t for t in InteriorTriangleIterator(dt)], fh) 

def test_small():
#     pts = [(3421275.7657, 3198467.4977), 
#            (3421172.5598, 3198197.546)
#            ]
#     triangulate(pts)
    # triangulate([(0,0), (0, -1)])
    #triangulate([(0,0)])
#     buggy = [(-120,90),(-60,40), (0,0),]# (-45, 35)]
#     triangulate(buggy)
#     triangulate([(0,0), (0,-20)])
    triangulate([(0,0), (10, 6), (6.5, 1.5), (2,3),(40,-20), (15, -4), (-120,90), (-60,40), (-45, 35)])
    #triangulate([(0,0), (-10, 6)])
#     triangulate([(0,0), (0, -6)])
