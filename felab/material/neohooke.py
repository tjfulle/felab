"""Neo Hookean HYPERELASTICITY, cannot be used for plane stress"""
import numpy as np


class neo_hooke(object):
    requires = ("nlgeom",)
    name = "Neo Hooke"

    def __init__(self, E, Nu):
        self.E, self.Nu = E, Nu

    def response(
        self,
        stress,
        statev,
        strain,
        dstrain,
        time,
        dtime,
        temp,
        dtemp,
        predef,
        dpred,
        ndir,
        nshr,
        ntens,
        coords,
        F0,
        F1,
        noel,
        kstep,
        kframe,
    ):

        if ndir != 3:
            raise TypeError(
                "Neo Hooke model requires " "three direct stress components"
            )

        # ELASTIC PROPERTIES
        C10 = self.E / (4.0 * (1.0 + self.Nu))
        D1 = 6.0 * (1.0 - 2.0 * self.Nu) / self.E

        # JACOBIAN AND DISTORTION TENSOR
        jac = F1[0, 0] * F1[1, 1] * F1[2, 2] - F1[0, 1] * F1[1, 0] * F1[2, 2]
        if nshr == 3:
            jac += F1[0, 1] * F1[1, 2] * F1[2, 0] + F1[0, 2] * F1[2, 1] * F1[1, 0]
            jac -= F1[0, 2] * F1[2, 0] * F1[1, 1] + F1[1, 2] * F1[2, 1] * F1[0, 0]
        scale = jac ** (-1.0 / 3.0)
        Fb = scale * F1

        # DEVIATORIC LEFT CAUCHY-GREEN DEFORMATION TENSOR
        Bb = np.zeros(ntens)
        Bb[0] = Fb[0, 0] * Fb[0, 0] + Fb[0, 1] * Fb[0, 1] + Fb[0, 2] * Fb[0, 2]
        Bb[1] = Fb[1, 0] * Fb[1, 0] + Fb[1, 1] * Fb[1, 1] + Fb[1, 2] * Fb[1, 2]
        Bb[2] = Fb[2, 0] * Fb[2, 0] + Fb[2, 1] * Fb[2, 1] + Fb[2, 2] * Fb[2, 2]
        Bb[3] = Fb[1, 0] * Fb[0, 0] + Fb[1, 1] * Fb[0, 1] + Fb[1, 2] * Fb[0, 2]
        if nshr == 3:
            Bb[4] = Fb[2, 0] * Fb[1, 0] + Fb[2, 1] * Fb[1, 1] + Fb[2, 2] * Fb[1, 2]
            Bb[5] = Fb[2, 0] * Fb[0, 0] + Fb[2, 1] * Fb[0, 1] + Fb[2, 2] * Fb[0, 2]

        I1B = sum(Bb[:3]) / 3.0
        G = 2.0 * C10 / jac
        K = 2.0 / D1 * (2.0 * jac - 1.0)
        Nu = 2.0 / D1 * (jac - 1.0)

        # CAUCHY STRESS
        stress[:ndir] = G * (Bb[:ndir] - I1B) + Nu
        stress[ndir : ndir + nshr] = G * Bb[ndir : ndir + nshr]

        # SPATIAL STIFFNESS
        D = np.zeros(ndir, ndir)
        G23 = G * 2.0 / 3.0
        D[0, 0] = G23 * (Bb[0] + I1B) + K
        D[0, 1] = -G23 * (Bb[0] + Bb[1] - I1B) + K
        D[0, 2] = -G23 * (Bb[0] + Bb[2] - I1B) + K
        D[0, 3] = G23 * Bb[3] / 2.0

        D[1, 1] = G23 * (Bb[1] + I1B) + K
        D[1, 2] = -G23 * (Bb[1] + Bb[2] - I1B) + K
        D[1, 3] = G23 * Bb[3] / 2.0

        D[2, 2] = G23 * (Bb[2] + I1B) + K
        D[2, 3] = -G23 * Bb[3]

        D[3, 3] = G * (Bb[0] + Bb[1]) / 2.0

        if nshr == 3:
            D[0, 4] = -G23 * Bb[4]
            D[0, 5] = G23 * Bb[5] / 2.0

            D[1, 4] = G23 * Bb[4] / 2.0
            D[1, 5] = -G23 * Bb[5]

            D[2, 4] = G23 * Bb[4] / 2.0
            D[2, 5] = G23 * Bb[5] / 2.0

            D[3, 4] = G * Bb[5] / 2.0
            D[3, 5] = G * Bb[4] / 2.0

            D[4, 4] = G * (Bb[0] + Bb[2]) / 2.0
            D[4, 5] = G * Bb[3] / 2.0

            D[5, 5] = G * (Bb[1] + Bb[2]) / 2.0

        return stress, D, statev
