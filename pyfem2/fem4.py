from numpy import *

from .constants import *
from .utilities import *
from .fem1 import FiniteElementModel
from .isoplib import IsoPElement

class Plane2DModel(FiniteElementModel):
    dimensions = 2

    def Solve(self, solver=None, **kwds):
        self.setup(IsoPElement)
        if solver is None:
            self.StaticPerturbation()
        elif solver == NEWTON:
            self.NewtonSolve(**kwds)
        else:
            raise ValueError('UNKNOWN SOLVER')

    def StaticPerturbation(self):
        z = zeros(self.numdof)
        K, rhs = self.assemble(self.dofs, z)
        Kbc, Fbc = self.apply_bc(K, rhs)
        self.dofs[:] = linsolve(Kbc, Fbc)
        R = dot(K, self.dofs) - rhs
        R = R.reshape(self.mesh.coord.shape)
        self.assemble(z, self.dofs, cflag=LP_OUTPUT)
        # CREATE NEW FRAME TO HOLD UPDATED STATE
        frame =self.steps.last.Frame(1.)
        self.advance(R=R)
        frame.converged = True

    def NewtonSolve(self, period=1., increments=5, maxiters=10, tolerance=1e-4,
                    relax=1.):

        time = array([0., self.steps.last.frames[-1].value])
        dtime = period / float(increments)

        istep = 1
        flags = [1, 0, 1, 0]
        for iframe in range(increments):

            load_fac = float(iframe+1) / float(increments)
            err1 = 1.

            # CREATE FRAME TO HOLD RESULTS
            frame = self.steps.last.Frame(dtime)

            u = zeros(self.numdof)

            # NEWTON-RAPHSON LOOP
            for nit in range(maxiters):

                K, rhs = self.assemble(self.dofs, u, time, dtime, istep, iframe+1,
                                       step_type=GENERAL, load_fac=load_fac)

                # Enforce displacement boundary conditions
                Kbc, Fbc = self.apply_bc(K, rhs, self.dofs, u, load_fac=load_fac)

                # --- Solve for the nodal displacement
                w = linsolve(Kbc, Fbc)

                # --- update displacement increment
                u += relax * w

                # --- Check convergence
                err1 = sqrt(dot(w, w) / dot(u, u))
                err2 = sqrt(dot(rhs, rhs)) / float(self.numdof)

                if err1 < tolerance or err2 < tolerance:
                    break

                # RESET DATA TO BEGINNING OF INCREMENT
                self.svars[1] = self.svars[0]

                continue

            else:
                raise RuntimeError('FAILED TO CONVERGE ON STEP {0}, '
                                   'FRAME {1}'.format(istep, iframe+1))

            time += dtime
            self.dofs += u
            self.advance()
            frame.converged = True

        return

    def advance(self, **kwds):
        self.svars[0] = self.svars[1]
        u = self.dofs.reshape(self.mesh.coord.shape)
        frame = self.steps.last.frames[-1]
        frame.field_outputs['U'].add_data(u)
        for (kwd, val) in kwds.items():
            frame.field_outputs[kwd].add_data(val)
        frame_n = self.steps.last.frames[frame.number-1]
        if not frame_n.converged:
            raise RuntimeError('ATTEMPTING TO UPDATE AN UNCONVERGED FRAME')
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            # GET STRAIN AND STRESS FROM LAST CONVERGED FRAME
            if not eb.eletyp.variables:
                continue
            ntens = eb.eletyp.ndir + eb.eletyp.nshr
            m = 1 if not eb.eletyp.integration else eb.eletyp.integration
            n = len(eb.eletyp.variables)
            for (e, xel) in enumerate(eb.labels):
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                ue = u[el.inodes]
                svars = self.svars[0,self.svtab[iel]].reshape(m,n,ntens)
                for (j, name) in enumerate(el.variables):
                    frame.field_outputs[eb.name,name].add_data(svars[:,j], e)
        return
