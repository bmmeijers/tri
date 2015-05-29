README
======

Tri - (Constrained) Delaunay Triangulation of Planar Straight Line Graphs


Installation
------------

- Download the source
- Run `python setup.py install` (or `python setup.py develop`)


Usage
-----

To triangulate a polygon with one outer and one inner shell::

    from tri import ToPointsAndSegments, triangulate
    from tri.delaunay import output_triangles, TriangleIterator

    # create points and segments for triangulation
    pts_segs = ToPointsAndSegments()
    pts_segs.add_polygon([[(0,0), (10,0), (5,10), (0,0)],
                          [(0,0), (8,2), (6,4), (5,7), (0,0)]
                          ],
                         )

    # triangulate the points and segments
    dt = triangulate(pts_segs.points, pts_segs.infos, pts_segs.segments)

    # write the output
    with open("tris.wkt", "w") as fh:
        output_triangles([t for t in TriangleIterator(dt)], fh)

The resulting file is readable with `QGIS <http://qgis.org>`_ (Add Delimited 
Text Layer).


Changelog
---------
See `CHANGES.txt`.


Bug reports
-----------
If you discover any bugs, feel free to create an issue on Bitbucket.

Please add as much information as possible to help us fixing the possible bug.
We also encourage you to help even more by forking and sending us a pull
request.

The issue tracker lives `here <https://bitbucket.org/bmmeijers/tri/issues>`_.


Maintainers
-----------

- `Martijn Meijers <https://bitbucket.org/bmmeijers>`_


Contributors
------------

- `Radan Šuba <https://bitbucket.org/radan/>`_


License
-------

`MIT License <https://www.tldrlegal.com/l/mit>`_