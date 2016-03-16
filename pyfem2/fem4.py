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
        self.advance(R=R)

    def NewtonSolve(self, period=1., increments=5, maxiters=10,
                    tolerance=1e-4, relax=1., tolerance1=1e-6):

        ti = self.steps.last.frames[-1].start
        time = array([ti, ti])
        dtime = period / float(increments)
        maxiter1 = max(int(maxiters/2.),1)

        istep = 1
        for iframe in range(increments):

            load_fac = float(iframe+1) / float(increments)
            err1 = 1.

            u = zeros(self.numdof)

            # NEWTON-RAPHSON LOOP
            for nit in range(maxiters):

                K, rhs = self.assemble(self.dofs, u, time, dtime, period,
                                       istep, iframe+1, ninc=nit+1,
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

                if nit < maxiter1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        logging.debug('CONVERGING TO LOSER TOLERANCE')
                        break

                continue

            else:
                raise RuntimeError('FAILED TO CONVERGE ON STEP {0}, '
                                   'FRAME {1}'.format(istep, iframe+1))

            time[1] += dtime
            self.dofs += u
            self.advance(dtime=dtime)

        return

    def advance(self, dtime=1., **kwds):

        frame_n = self.steps.last.frames[-1]
        if not frame_n.converged:
            raise RuntimeError('ATTEMPTING TO UPDATE AN UNCONVERGED FRAME')

        # ADVANCE STATE VARIABLES
        self.svars[0] = self.svars[1]
        u = self.dofs.reshape(self.mesh.coord.shape)

        # CREATE FRAME TO HOLD RESULTS
        frame = self.steps.last.Frame(dtime)
        frame.field_outputs['U'].add_data(u)
        for (kwd, val) in kwds.items():
            frame.field_outputs[kwd].add_data(val)

        for (ieb, eb) in enumerate(self.mesh.eleblx):
            if not eb.eletyp.variables:
                continue

            # PASS VALUES FROM SVARS TO THE FRAME OUTPUT
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

        frame.converged = True
        return
