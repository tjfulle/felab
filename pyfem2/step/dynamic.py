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
    def run(self):

        period, increments = self.period, self.increments

        time = array([0., self.start])
        dtime = period / float(increments)

        shape = (self.model.numdof, 2)
        un = zeros(shape)
        vn = zeros(shape)
        an = zeros(shape)

        count = 0

        m = mass_matrix
        F = external_force_array_at_time_equal_0
        an[:,0] = F / m

        a, b = self.alpha, self.beta
        for n in range(increments):

            dltyp, dload = self.dload(time[1])
            mass, rhs = self.model.assemble(self.dofs, u, Q,
                                            self.svtab, self.svars,
                                            dltyp, dload, self.predef,
                                            self.procedure, MASS_AND_RHS, time=time,
                                            istep=self.number, iframe=iframe+1)

            an[:,1] = rhs / mass
            vn[:,1] = vn[:,0] + ((1. - a) * an[:,0] + a * an[:,1]) * dtime
            un[:,1] = un[:,0] + (vn[:,0]+.5*((1.-b)*an[:,0]+b*an[:,1])*dtime)*dtime

            un[:,0] = un[:,1]
            vn[:,0] = vn[:,1]
            an[:,0] = an[:,1]

            time[1] += dtime

            self.dofs[:] = un[:,0]
            self.vel[:] = vn[:,0]
            self.accel[:] = an[:,0]
            self.advance(dtime, self.dofs)

        return
