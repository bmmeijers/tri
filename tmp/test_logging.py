# myapp.py
import logging
from tri import ToPointsAndSegments, triangulate
from tri.delaunay import TriangleIterator



def main():
    logging.basicConfig(filename='test_logging.log', level=logging.DEBUG)
    
    pts_segs = ToPointsAndSegments()
    pts_segs.add_polygon([[(0,0), (10,0), (5,10), (0,0)],
                          [(0,0), (8,2), (6,4), (5,7), (0,0)]
                          ],
                         )
    dt = triangulate(pts_segs.points, pts_segs.infos, pts_segs.segments)
    

if __name__ == '__main__':
    main()