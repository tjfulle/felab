from numpy import *

from ..constants import *
from ..utilities import *
from ._step import SDStep

class DynamicStep(SDStep):
    procedure = DYNAMIC
    def __init__(self, model, number, name, previous, period, increments, nlgeom):
        super(DynamicStep, self).__init__(model, number, name, previous, period)
        self.nlgeom = nlgeom
        self.increments = increments

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, alpha=.5, beta=0.):

        period, increments = self.period, self.increments

        time = array([0., self.start])
        dtime = period / float(increments)

        shape = (self.model.numdof, 2)
        un = zeros(shape)
        vn = zeros(shape)
        an = zeros(shape)

        # GET MASS AND RHS AT TIME 0
        Q = self.cload(time[0])
        dltyp, dload = self.dload(time[0])

        mass, rhs = self.model.assemble(
            self.dofs, un[:,0], Q, self.svtab, self.svars, dltyp, dload,
            self.predef, self.procedure, GENERAL, cflag=MASS_AND_RHS, time=time)

        an[:,0] = rhs / mass

        for n in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            mass, rhs = self.model.assemble(
                self.dofs, un[:,0], Q, self.svtab, self.svars, dltyp, dload,
                self.predef, self.procedure, GENERAL, cflag=MASS_AND_RHS,
                time=time, istep=self.number, iframe=n+1)

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
