'''
Created on Nov 13, 2018

@author: martijn
'''


def output_vertices_iter_info(V, fh):
    """Output list of vertices as WKT to text file (for QGIS)"""
    fh.write("id;wkt;finite;{}\n".format(";".join(
        "info_{}".format(i) for i, _ in enumerate(V[0].info, start=1))))
    for v in V:
        fh.write("{0};POINT({1});{2};{3}\n".format(
            id(v), v, v.is_finite, ";".join(str(each) for each in v.info)))


def output_vertices(V, fh):
    """Output list of vertices as WKT to text file (for QGIS)"""
    fh.write("id;wkt;finite;info\n")
    for v in V:
        fh.write("{0};POINT({1});{2};{3}\n".format(
            id(v), v, v.is_finite, v.info))


def output_triangles(T, fh):
    """Output list of triangles as WKT to text file (for QGIS)"""
    fh.write("id;wkt;n0;n1;n2;v0;v1;v2;c0;c1;c2;finite;info\n")
    for t in T:
        if t is None:
            continue
        fh.write("{0};{1};"
                 "{2[0]};{2[1]};{2[2]};"
                 "{3[0]};{3[1]};{3[2]};"
                 "{4[0]};{4[1]};{4[2]};"
                 "{5};{6}\n".format(
                    id(t), t,
                    [id(n) for n in t.neighbours],
                    [id(v) for v in t.vertices],
                    t.constrained,
                    t.is_finite, t.info))


def output_edges(E, fh):
    fh.write("id;side;wkt\n")
    for e in E:
        fh.write("{0};{1};"
                 "LINESTRING({2[0][0]} {2[0][1]}, {2[1][0]} {2[1][1]})".format(
                    id(e.triangle), e.side,
                    e.segment))
        fh.write("\n")
