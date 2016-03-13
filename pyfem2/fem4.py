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
            self.NewtonSolve()
        else:
            raise ValueError('UNKNOWN SOLVER')

    def StaticPerturbation(self):
        flags = [1, 0, 1, 1]
        du = zeros(self.numdof)
        K, rhs = self.assemble(self.dofs, du, [0, 0], 1., 1, 1, flags)
        Kbc, Fbc = self.apply_bc(K, rhs)
        self.dofs[:] = linsolve(Kbc, Fbc)
        R = dot(K, self.dofs) - rhs
        R = R.reshape(self.mesh.coord.shape)
        u = self.dofs.reshape(self.mesh.coord.shape)

        # Create new frame to hold updated state
        frame =self.steps.last.Frame(1.)
        self.update_state(u, 1., R=R)
        frame.converged = True

    def NewtonSolve(self, period=1., increments=5, maxiters=10, tolerance=1e-4,
                    relax=1.):

        time = self.steps.last.frames[-1].value
        dtime = period / float(increments)

        istep = 1
        for iframe in range(1, increments+1):
            load_fac = float(iframe) / float(increments)
            err1 = 1.

            # create frame to hold results
            frame = self.steps.last.Frame(dtime)

            u = zeros(self.numdof)

            # Newton-Raphson loop
            for nit in range(maxiters):

                K, F, Q = self.assemble(self.dofs, u, time, dtime, istep, iframe)
                R = self.assemble_global_residual(u, dtime)

                # Compute global force
                rhs = load_fac * (F + Q) - R

                # Enforce displacement boundary conditions
                Kbc, Fbc = self.apply_bc(K, rhs, self.dofs, u, load_fac)

                # --- Solve for the nodal displacement
                du = linsolve(Kbc, rhs)

                # --- update displacement increment
                u += relax * du

                # --- Check convergence
                dusq = dot(u, u)
                err1 = dot(du, du)
                if dusq != 0.:
                    err1 = sqrt(err1 / dusq)
                err2 = sqrt(dot(rhs, rhs)) / float(self.numdof)

                self.update_state(u, dtime)

                if err1 < tolerance or err2 < tolerance:
                    break

                continue

            else:
                raise RuntimeError('FAILED TO CONVERGE ON STEP {0}, '
                                   'FRAME {1}'.format(istep, iframe))

            frame.converged = True
            time += dtime
            self.dofs += u

        return

    def update_state(self, u, dtime, **kwds):
        frame = self.steps.last.frames[-1]
        frame.field_outputs['U'].add_data(u)
        for (kwd, val) in kwds.items():
            frame.field_outputs[kwd].add_data(val)
        frame_n = self.steps.last.frames[frame.number-1]
        if not frame_n.converged:
            raise RuntimeError('ATTEMPTING TO UPDATE AFTER UNCONVERGED FRAME')
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            # get strain and stress from last converged frame
            E = frame_n.field_outputs[eb.name, 'E']
            S = frame_n.field_outputs[eb.name, 'S']
            for (e, xel) in enumerate(eb.labels):
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                ue = u[el.inodes]
                de1, e1, s1 = el.update_state(ue, E.data[e], S.data[e])
                frame.field_outputs[eb.name, 'S'].add_data(s1, e)
                frame.field_outputs[eb.name, 'E'].add_data(e1, e)
                frame.field_outputs[eb.name, 'DE'].add_data(de1, e)
        return
