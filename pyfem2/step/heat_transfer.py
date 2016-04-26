from numpy import *

from ..constants import *
from ..utilities import *
from ..elemlib import PlaneDiffussiveHeatTransferTria3
from ._step import Step

class HeatTransferStep(Step):
    procedure = HEAT_TRANSFER
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
            self.assign_hsrc(ielem, a[ix])

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self):
        frame = self.Frame()
        time = array([0., self.start])
        du = zeros(self.model.numdof)
        qe = zeros_like(self.dofs)
        dltyp, dload = self.dload(self.period)
        X = self.dofvals(self.period)
        K, rhs = self.model.assemble(self.dofs, du, qe, frame.field_outputs,
                                     dltyp, dload, self.predef,
                                     self.procedure, DIRECT, time=time)
        Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, X)
        self.dofs[:] = linsolve(Kbc, Fbc)
        react = dot(K, self.dofs) - rhs
        u, R, temp = self.model.format_dof(self.dofs)
        RF, M, Q = self.model.format_dof(react)
        frame.snapshot(self.period, T=temp, Q=Q)
