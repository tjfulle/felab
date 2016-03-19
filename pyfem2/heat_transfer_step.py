from numpy import *

from .constants import *
from .utilities import *
from .step import Step
from .elemlib import PlaneDiffussiveHeatTransferTria3

class HeatTransferStep(Step):
    def __init__(self, model, number, name, previous, period):
        super(HeatTransferStep, self).__init__(model, number, name, previous,
                                               period)

        # CHECK ELEMENTS
        eletyp = (PlaneDiffussiveHeatTransferTria3,)
        if not all([isinstance(el, eletyp) for el in self.model.elements]):
            raise UserInputError('INCORRECT ELEMENT TYPE FOR HEAT TRANSFER STEP')

    # ----------------------------------------------------------------------- #
    # --- HEAT TRANSFER LOADINGS -------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def SurfaceFlux(self, surface, qn):
        surf = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.dltyp[iel].append(SFLUX)
            self.dload[iel].append([iedge, qn])

    def SurfaceConvection(self, surface, Too, h):
        if self.model.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        surf = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.dltyp[iel].append(SFILM)
            self.dload[iel].append([iedge, Too, h])

    def HeatGeneration(self, region, amplitude):
        if region == ALL:
            xelems = sorted(self.model.mesh.elemap.keys())
        else:
            xelems = region
        inodes = []
        for eb in self.model.mesh.eleblx:
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
            self.dltyp[ielem].append(HSRC)
            self.dload[ielem].append(a[ix])

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self):
        du = zeros(self.model.numdof)
        qe = zeros_like(self.dofs)
        K, rhs = self.model.assemble(self.dofs, du, qe, self.svtab, self.svars,
                                     self.dltyp, self.dload, self.predef,
                                     HEAT_TRANSFER, DIRECT)
        Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, self.dofvals)
        self.dofs[:] = linsolve(Kbc, Fbc)
        react = dot(K, self.dofs) - rhs
        self.advance(self.period, self.dofs, react)

        # CREATE NEW FRAME TO HOLD UPDATED STATE
        #frame = self.Frame(self.period)
        #frame.field_outputs['T'].add_data(self.dofs)
        #frame.field_outputs['Q'].add_data(Q)
        #frame.converged = True
        #self.T = self.dofs
