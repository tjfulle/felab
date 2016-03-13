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
        K, F, Q = self.assemble()
        Kbc, Fbc = self.apply_bc(K, F+Q)
        u = linsolve(Kbc, Fbc)
        Ft = dot(K, u)
        self.dofs[:] = u.reshape(self.dofs.shape)
        R = Ft - F - Q
        R = R.reshape(self.mesh.coord.shape)
        u = u.reshape(self.mesh.coord.shape)

        # Create new frame to hold updated state
        frame =self.steps.last.Frame(1.)
        self.update_state(u, 1., R=R)
        frame.converged = True

    def NewtonSolve(self, period=1., increments=5, maxiters=10, tolerance=1e-4,
                    relax=1.):

        u = zeros(self.numdof)

        time = self.steps.last.frames[-1].value
        dtime = period / float(increments)

        istep = 1
        for iframe in range(1, increments+1):
            load_fac = float(iframe) / float(increments)
            err1 = 1.

            # create frame to hold results
            frame = self.steps.last.Frame(dtime)

            # Newton-Raphson loop
            for nit in range(maxiters):

                K, F = self.assemble(u, du, time, dtime, istep, iframe)

                K = self.assemble_global_stiffness(u, dtime)
                F, Q = self.assemble_global_force(self.dload, self.sload)
                R = self.assemble_global_residual(u, dtime)

                # Compute global force
                rhs = load_fac * (F + Q) - R

                # Enforce displacement boundary conditions
                Kbc, Fbc = self.apply_bc(K, rhs, u, load_fac)

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
