from numpy import *

from ..utilities import *
from .element import Element

# --------------------------------------------------------------------------- #
# ------------------------ HEAT TRANSFER ELEMENT ---------------------------- #
# --------------------------------------------------------------------------- #
class CHTElement(Element):
    edges = []
    nodes = None
    gaussp = None
    elefab = None
    gaussw = None
    signature = None
    dimensions = None

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

    def response(self, u, du, time, dtime, istep, iframe, svars, dltyp, dload,
                 predef, procedure, nlgeom, cflag, step_type, load_fac):

        # --- ELEMENT STIFFNESS AND FORCE

        compute_stiff = cflag in (STIFF_AND_FORCE, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_FORCE, FORCE_ONLY)

        if compute_stiff:
            Ke = self.conduction_stiff_contrib()

        if compute_force:
            Fe = zeros(3)

        for (i, typ) in enumerate(dltyp):
            # CONTRIBUTIONS FROM EXTERNAL LOADS

            if typ == SFILM:
                # EVALUATE THE CONVECTION CONTRIBUTION
                iedge, Too, h = dload[i]
                if compute_stiff:
                    Ke += self.convection_stiff_contrib(iedge, h)
                if compute_force:
                    Fe += self.convection_flux_array(iedge, Too, h)

            elif typ == HSRC:
                # EVALUATE THE HEAT SOURCE CONTRIBUTION
                if compute_force:
                    Fe += self.heat_source(dload[i])

            elif typ == SFLUX:
                # EVALUATE THE BOUNDARY FLUX CONTRIBUTION
                if compute_force:
                    iedge, qn = dload[i]
                    Fe += self.conduction_flux_array(iedge, qn)

        if cflag == STIFF_AND_FORCE:
            return Ke, Fe

        elif cflag == STIFF_ONLY:
            return Ke

        elif cflag == FORCE_ONLY:
            return Fe
