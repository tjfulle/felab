from numpy import *

from .constants import *
from .utilities import *
from .fem1 import FiniteElementModel
from .isoplib import IsoPElement

class Plane2DModel(FiniteElementModel):
    numdim = 2

    def Solve(self, solver=None, **kwds):
        self.setup(IsoPElement)
        if solver is None:
            self.StaticPerturbation()
        elif solver == NEWTON:
            self.NewtonSolve()
        else:
            raise ValueError('Unknown solver')

    def StaticPerturbation(self):
        K = self.assemble_global_stiffness()
        F, Q = self.assemble_global_force(self.dload, self.sload)
        Kbc, Fbc = self.apply_bc(K, F+Q)
        u = linsolve(Kbc, Fbc)
        Ft = dot(K, u)
        self.dofs[:] = u.reshape(self.dofs.shape)
        R = Ft - F - Q
        R = R.reshape(self.mesh.coord.shape)
        u = u.reshape(self.mesh.coord.shape)

        # Create new frame to hold updated state
        frame =self.steps.last.Frame(1.)
        self.update_state(u, 1.)
        frame.field_outputs['R'].add_data(R)
        frame.converged = True

    def NewtonSolve(self, period=1., increments=5, maxiters=10, tolerance=1e-4,
                    relax=1.):

        u = zeros(self.numdof)
        w = zeros(self.numdof)

        time = self.steps.last.frames[-1].value
        dtime = period / float(increments)

        istep = 1
        for iframe in range(1, increments+1):
            load_fac = float(iframe) / float(increments)
            err1 = 1.

            # create frame to hold results
            self.steps.last.Frame(dtime)

            # Newton-Raphson loop
            for nit in range(maxiters):

                K = self.assemble_global_stiffness(w, dtime)
                F, Q = self.assemble_global_force(self.dload, self.sload)
                R = self.assemble_global_residual(w, dtime)

                # Compute global force
                rhs = load_fac * (F + Q) - R

                # Enforce displacement boundary conditions
                Kbc, Fbc = self.apply_bc(K, rhs, w, load_fac)

                # --- Solve for the nodal displacement
                dw = linsolve(Kbc, rhs)

                # --- update displacement increment
                w += relax * dw

                # --- Check convergence
                dusq = dot(w, w)
                err1 = dot(dw, dw)
                if dusq != 0.:
                    err1 = sqrt(err1 / dusq)
                err2 = sqrt(dot(rhs, rhs)) / float(self.numdof)

                if err1 < tolerance or err2 < tolerance:
                    break

                continue

            else:
                raise RuntimeError('Failed to converge on step {0}, '
                                   'frame {1}'.format(istep, iframe))

            time += dtime
            u += w
            self.update_state(u, dtime)

        return

    def update_state(self, u, dtime):
        frame = self.steps.last.frames[-1]
        frame.field_outputs['U'].add_data(u)
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            E = frame.field_outputs[eb.name, 'E']
            S = frame.field_outputs[eb.name, 'S']
            for (e, xel) in enumerate(eb.labels):
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                ue = u[el.nodes]
                de1, e1, s1 = el.update_state(ue, E.data[e], S.data[e])
                frame.field_outputs[eb.name, 'S'].add_data(s1,e)
                frame.field_outputs[eb.name, 'E'].add_data(e1,e)
                frame.field_outputs[eb.name, 'DE'].add_data(de1,e)

        return

    def update_kinematics(self, u, dtime):
        """Update kinematics to dtime"""
        frame = self.steps.last.frames[-1]
        frame.field_outputs['U'].add_data(u)
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            E = frame.field_outputs[eb.name, 'E']
            for (e, xel) in enumerate(eb.labels):
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                de1, e1 = el.update_kinematics(ue, dtime, E)
                frame.field_outputs[eb.name, 'E'].add_data(e1, e)
                frame.field_outputs[eb.name, 'DE'].add_data(de1, e)

        return
