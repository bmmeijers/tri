from tri.delaunay import ToPointsAndSegments, triangulate


from tri.delaunay.insert_kd import kdsort
import unittest

#class TestTri(unittest.TestCase):
#    def test_convert(self):
#        conv = ToPointsAndSegments()
#        conv.add_polygon(
#            [[(0, 0), (22, 0), (14, 10), (2, 8), (0, 6.5), (0, 0)]], info = 'a')
#        assert len(conv.points) == 5
#        assert len(conv.segments) == 5
#        assert len(conv.infos) == 5
#        # all infos for the vertices should be a list with 'a'
#        for item in conv.infos:
#            assert item[1] == ['a']

#    def test_triangulate(self):
#        conv = ToPointsAndSegments()
#        conv.add_polygon(
#            [[(0, 0), (22, 0), (14, 10), (2, 8), (0, 6.5), (0, 0)]], info = 'a')
#        dt = triangulate(conv.points)
#        self.assertEqual(len(dt.triangles), 11)
#        self.assertEqual(len(dt.vertices), len(conv.points))
#        self.assertEqual(len(dt.vertices), 5)

#    def test_triangulate_with_infos(self):
#        conv = ToPointsAndSegments()
#        conv.add_polygon(
#            [[(0, 0), (22, 0), (14, 10), (2, 8), (0, 6.5), (0, 0)]], info = 'a')
#        dt = triangulate(conv.points, conv.infos)
#        self.assertEqual(len(dt.triangles), 11)
#        self.assertEqual(len(dt.vertices), len(conv.points))
#        self.assertEqual(len(dt.vertices), 5)
#        
#    def test_triangulate_with_infos_and_segments(self):
#        conv = ToPointsAndSegments()
#        conv.add_polygon(
#            [[(0, 0), (22, 0), (14, 10), (2, 8), (0, 6.5), (0, 0)]], info = 'a')
#        dt = triangulate(conv.points, conv.infos, conv.segments)
#        self.assertEqual(len(dt.triangles), 11)
#        self.assertEqual(len(dt.vertices), len(conv.points))
#        self.assertEqual(len(dt.vertices), 5)
#        for v, i in zip(dt.vertices, [3, 4, 0, 1, 2]):
#            assert(v.info[0] == i)
#            assert(v.info[1] == ['a'])

class TestTriSimple(unittest.TestCase):
#    def test_triangle(self):
#        conv = ToPointsAndSegments()
#        conv.add_point((10, 0), info = 'a')
#        conv.add_point((-2, 8), info = 'b')
#        conv.add_point((-2, -8) , info = 'c')
#        conv.add_segment((10, 0), (-2, 8))
#        conv.add_segment((-2, 8), (-2, -8))
#        conv.add_segment((-2, -8), (10, 0))
#        dt = triangulate(conv.points, conv.infos, conv.segments)
#        self.assertEqual(len(dt.triangles), 7)
#        self.assertEqual(len(dt.vertices), 3)


    def test_dent(self):
#        conv = ToPointsAndSegments()
#        polygon = [[(0, 0), (10., 0), (10,20), (-0.5,20.), (-0.5,11.), (-1,11), (-1,10), (0,10), (0,0)]]
#        conv.add_polygon(polygon)
#        print('pts', conv.points)
#        print('segments', conv.segments)

#        print('pts - sorted', kdsort(conv.points))
#        dt = triangulate(conv.points, conv.infos, conv.segments, output = True)

        # ccw around shape
        points = [(0.0, 0.0), (10.0, 0.0), (10.0, 20.0), (-0.5, 20.0), (-0.5, 11.0), (-1.0, 11.0), (-1.0, 10.0), (0.0, 10.0)]
        infos = ['a','b','c','d',  'e' ,'f','g','h']
        segments = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)]
        dt = triangulate(points, infos, segments, output = True)


if __name__ == "__main__":
    import logging
    logging.getLogger().setLevel(logging.DEBUG)
    unittest.main()    
