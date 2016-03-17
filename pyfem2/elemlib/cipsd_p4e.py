from numpy import zeros

from ._cipsdq4 import CIPDSDQ4Element

# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# --------------------------------------------------------------------------- #
class PlaneStrainQuad4(CIPDSDQ4Element):
    ndir, nshr = 3, 1
    def bmatrix(self, dN):
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        return B
