'''
Created on Nov 13, 2018

@author: martijn
'''

import warnings

try:
    from geompreds import orient2d, incircle
except ImportError:
    warnings.warn(
        "Robust predicates not available, falling back on non-robust ones."
    )

    def orient2d(pa, pb, pc):
        """Direction from pa to pc, via pb, where returned value is as follows:

        left:     + [ = ccw ]
        straight: 0.
        right:    - [ = cw ]

        returns twice signed area under triangle pa, pb, pc
        """
        detleft = (pa[0] - pc[0]) * (pb[1] - pc[1])
        detright = (pa[1] - pc[1]) * (pb[0] - pc[0])
        det = detleft - detright
        return det

    def incircle(pa, pb, pc, pd):
        """Tests whether pd is in circle defined by the 3 points pa, pb and pc
        """
        adx = pa[0] - pd[0]
        bdx = pb[0] - pd[0]
        cdx = pc[0] - pd[0]
        ady = pa[1] - pd[1]
        bdy = pb[1] - pd[1]
        cdy = pc[1] - pd[1]
        bdxcdy = bdx * cdy
        cdxbdy = cdx * bdy
        alift = adx * adx + ady * ady
        cdxady = cdx * ady
        adxcdy = adx * cdy
        blift = bdx * bdx + bdy * bdy
        adxbdy = adx * bdy
        bdxady = bdx * ady
        clift = cdx * cdx + cdy * cdy
        det = alift * (bdxcdy - cdxbdy) + \
            blift * (cdxady - adxcdy) + \
            clift * (adxbdy - bdxady)
        return det
