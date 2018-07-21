from numpy import *
from numpy.linalg import inv

class linear_elastic(object):
    """Linear elastic material """
    name = 'elastic'
    def __init__(self, Lambda, Mu):
        self.Lambda, self.Mu = Lambda, Mu
    def response(self, stress, statev, strain, dstrain, time, dtime,
                 temp, dtemp, predef, dpred, ndir, nshr, ntens,
                 coords, F0, F, noel, kstep, kinc):

        C11 = self.Lambda + 2*self.Mu
        C12 = self.Lambda
        C44 = self.Mu
        D = array([[C11, C12, C12, 0,   0,   0  ],
                   [C12, C11, C12, 0,   0,   0  ],
                   [C12, C12, C11, 0,   0,   0  ],
                   [0,   0,   0,   C44, 0,   0  ],
                   [0,   0,   0,   0,   C44, 0  ],
                   [0,   0,   0,   0,   0,   C44]])

        if nshr == 1:
            # MODIFY THE STIFFNESS FOR 2D ACCORDING TO:
            # 1) PLANE STRAIN: REMOVE ROWS AND COLUMNS OF THE STIFFNESS
            #    CORRESPONDING TO THE PLANE OF ZERO STRAIN
            # 2) PLANE STRESS: INVERT THE STIFFNESS AND REMOVE THE ROWS
            #    AND COLUMNS OF THE COMPLIANCE CORRESPONDING THE PLANE OF
            #    ZERO STRESS.
            if ndir == 2:
                # PLANE STRESS
                # INVERT THE STIFFNESS TO GET THE COMPLIANCE
                idx = [[[0], [1], [3]], [0, 1, 3]]
                D = inv(inv(D)[idx])
            elif ndir == 3:
                # PLANE STRAIN
                idx = [[[0], [1], [2], [3]], [0, 1, 2, 3]]
                D = D[idx]

        stress += dot(D, dstrain)
        return stress, statev, D
