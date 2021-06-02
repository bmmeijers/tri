"""Tri - Constrained Delaunay Triangulation of Planar Straight Line Graphs
"""

import logging
import time
from datetime import datetime

from tri.delaunay.insert_kd import triangulate
from tri.delaunay.helpers import ToPointsAndSegments
# from tri.delaunay.insert_hcpo import triangulate


__version__ = '0.3.1.dev0'
__license__ = 'MIT License'
__author__ = 'Martijn Meijers'
__all__ = ("triangulate", "ToPointsAndSegments")


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    from tri.delaunay.helpers import random_circle_vertices
    pts = random_circle_vertices(150000)
    triangulate(pts)
