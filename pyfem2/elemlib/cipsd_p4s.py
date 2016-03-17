from numpy import zeros

from ._cipsdq4 import CIPDSDQ4Element

# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRESS ELEMENT ----------------------- #
# --------------------------------------------------------------------------- #
class PlaneStressQuad4(CIPDSDQ4Element):
    ndir, nshr = 2, 1
    def bmatrix(self, dN):
        B = zeros((3, 8))
        B[0, 0::2] = B[2, 1::2] = dN[0, :]
        B[1, 1::2] = B[2, 0::2] = dN[1, :]
        return B
