import logging
from numpy import *
from numpy.linalg import det, inv

from ..constants import *
from ..utilities import *
from .isop_base import isop_base


class DCMDN(isop_base):

    @staticmethod
    def variables():
        return None

    def conduction_stiff_contrib(self):
        raise NotImplementedError

    def convection_stiff_contrib(self, edge, h):
        raise NotImplementedError

    def heat_source(self, f):
        raise NotImplementedError

    def conduction_flux_array(self, edge, qn):
        raise NotImplementedError

    def convection_flux_array(self, edge, Too, h):
        raise NotImplementedError

    def boundary_flux_array(self, edge, qn):
        raise NotImplementedError

    def response(self, rhs, A, svars, energy, u, du, v, a, time, dtime,
                 kstage, kinc, dltyp, dlmag, predef, lflags,
                 ddlmag, mdload, pnewdt):

        if lflags[2] not in (STIFF_AND_RHS, STIFF_ONLY, RHS_ONLY):
            raise NotImplementedError

        # --- ELEMENT STIFFNESS AND FORCE
        if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY):
            A[:] = self.conduction_stiff_contrib()

        if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
            rhs[:] = 0.

        for (i, typ) in enumerate(dltyp):
            # CONTRIBUTIONS FROM EXTERNAL LOADS

            if typ == SFILM:
                # EVALUATE THE CONVECTION CONTRIBUTION
                iedge, Too, h = dlmag[i]
                if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY):
                    A[:] += self.convection_stiff_contrib(iedge, h)
                if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
                    rhs[:] += self.convection_flux_array(iedge, Too, h)

            elif typ == HSRC:
                # EVALUATE THE HEAT SOURCE CONTRIBUTION
                if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
                    rhs[:] += self.heat_source(dlmag[i])

            elif typ == SFLUX:
                # EVALUATE THE BOUNDARY FLUX CONTRIBUTION
                if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
                    iedge, qn = dlmag[i]
                    rhs[:] += self.conduction_flux_array(iedge, qn)
