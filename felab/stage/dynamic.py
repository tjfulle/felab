from numpy import *

from .stage import sd_stage
from ..constants import *
from ..utilities import *
from ..assembly import assemble_system

class dynamic_stage(sd_stage):
    def __init__(self, model, number, name, previous, period, **kwds):
        super(dynamic_stage, self).__init__(model, number, name, previous, period)
        for (key, val) in kwds.items():
            if key == 'increments':
                key = '_increments'
            setattr(self, key, val)

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, period=1., increments=10, alpha=.5, beta=0.):

        period = getattr(self, 'period', period)
        increments = getattr(self, '_increments', increments)

        time = array([0., self.start])
        dtime = period / float(increments)

        shape = (self.model.numdof, 2)
        un = zeros(shape)
        vn = zeros(shape)
        an = zeros(shape)

        # GET MASS AND RHS AT TIME 0
        Q = self.cload(time[0])
        dltyp, dlmag = self.dload(time[0])

        energy = v = a = ddlmag = mdload = None
        pnewdt = zeros(1)
        lflags = [DYNAMIC, SMALL_DISPLACEMENT, MASS_AND_RHS, GENERAL, 0]
        A = zeros((self.model.numdof, self.model.numdof))
        rhs = zeros(self.model.numdof)
        assemble_system(rhs, A, self.svtab, self.svars, energy, Q,
                        self.dofs, un[:,0], vn[:,0], an[:,0], time, dtime,
                        0, 1, 1, dltyp, dlmag, self.predef, lflags,
                        ddlmag, mdload, pnewdt, period,
                        self.model.mesh.element_blocks,
                        self.model.mesh.elemap, self.model.elements,
                        self.model.eftab)

        # LUMPED MASS MATRIX
        mass = array([sum(row) for row in A])
        an[:,0] = rhs / mass

        for kinc in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            rhs = zeros(self.model.numdof)
            A = zeros((self.model.numdof, self.model.numdof))
            assemble_system(rhs, A, self.svtab, self.svars, energy, Q,
                            self.dofs, un[:,0], vn[:,0], an[:,0], time, dtime,
                            self.number, kinc+1, 1, dltyp, dlmag, self.predef,
                            lflags, ddlmag, mdload, pnewdt, self.period,
                            self.model.mesh.element_blocks,
                            self.model.mesh.elemap, self.model.elements,
                            self.model.eftab)
            # LUMPED MASS MATRIX
            mass = array([sum(row) for row in A])

            # UPDATE ACCELERATION
            an[:,1] = rhs / mass

            # UPDATE VELOCITY
            dv = ((1. - alpha) * an[:,0] + alpha * an[:,1]) * dtime
            vn[:,1] = vn[:,0] + dv

            # UPDATE DISPLACEMENT
            du = (vn[:,0] + .5*((1.-beta)*an[:,0] + beta*an[:,1])*dtime)*dtime
            un[:,1] = un[:,0]  + du
            un[self.doftags,1] = X

            # ADVANCE KINEMATICS
            un[:,0] = un[:,1]
            vn[:,0] = vn[:,1]
            an[:,0] = an[:,1]

            time += dtime

            self.dofs[:] = un[:,0]
            #self.vel[:] = vn[:,0]
            #self.accel[:] = an[:,0]

            self.advance(dtime, self.dofs)

        return
