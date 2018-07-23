from numpy import *

from .stage import load_stage
from ..x.constants import *
from ..x.utilities import *
from ..elemlib import DC2D3
from ..assembly import assemble_system, apply_boundary_conditions

class diffusive_ht_stage(load_stage):
    def __init__(self, model, number, name, previous, period):
        super(diffusive_ht_stage, self).__init__(model, number, name, previous, period)

        # CHECK ELEMENTS
        eletyp = (DC2D3,)
        if not all([isinstance(el, eletyp) for el in self.model.elements]):
            raise UserInputError('INCORRECT ELEMENT TYPE FOR HEAT TRANSFER STAGE')

    # ----------------------------------------------------------------------- #
    # --- HEAT TRANSFER LOADINGS -------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def SurfaceFlux(self, surface, qn):
        surf = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.assign_sflux(iel, iedge, qn)

    def SurfaceConvection(self, surface, Too, h):
        if self.model.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        surf = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.assign_sfilm(iel, iedge, Too, h)

    def HeatGeneration(self, region, amplitude):
        if region == ALL:
            xelems = sorted(self.model.mesh.elemap.keys())
        else:
            xelems = region
        inodes = []
        for eb in self.model.mesh.element_blocks:
            for (i, xel) in enumerate(eb.labels):
                if xel in xelems:
                    inodes.extend(eb.elecon[i])
        inodes = unique(inodes)
        if hasattr(amplitude, '__call__'):
            # AMPLITUDE IS A FUNCTION
            x = self.model.mesh.coord[inodes]
            a = amplitude(x)
        elif not is_listlike(amplitude):
            a = amplitude * ones(len(inodes))
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError('HEAT GENERATION AMPLITUDE MUST HAVE '
                                     'LENGTH {0}'.format(len(inodes)))
            a = asarray(amplitude)
        nodmap = dict(zip(inodes, range(inodes.shape[0])))
        for xelem in xelems:
            ielem = self.model.mesh.elemap[xelem]
            ix = [nodmap[n] for n in self.model.elements[ielem].inodes]
            self.assign_hsrc(ielem, a[ix])

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self):
        # define arguments needed for assembly
        time = array([0., self.start])
        du = zeros(self.model.numdof)
        Q = zeros_like(self.dofs)
        dltyp, dlmag = self.dload(self.period)
        X = self.dofvals(self.period)
        energy = v = a = ddlmag = mdload = pnewdt = None
        time = array([1., 1.])
        dtime = 1.
        lflags = [HEAT_TRANSFER_STEADY_STATE, SMALL_DISPLACEMENT, STIFF_AND_RHS, GENERAL, 0]

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        rhs = zeros(self.model.numdof)
        K = zeros((self.model.numdof, self.model.numdof))
        assemble_system(rhs, K, self.svtab, self.svars, energy, Q,
                        self.dofs, du, v, a, time, dtime, self.number, 1, 0,
                        dltyp, dlmag, self.predef, lflags,
                        ddlmag, mdload, pnewdt, self.period,
                        self.model.mesh.element_blocks, self.model.mesh.elemap,
                        self.model.elements, self.model.eftab)
        self._K = K
        Kbc, Fbc = apply_boundary_conditions(K, rhs, self.doftags, X)
        self.dofs[:] = linsolve(Kbc, Fbc)
        react = dot(K, self.dofs) - rhs
        self.advance(self.period, self.dofs, react)
