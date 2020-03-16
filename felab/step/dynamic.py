import numpy as np

from felab.step.step import sd_step
from felab.assembly import assemble_system
from felab.constants import DYNAMIC, SMALL_DISPLACEMENT, MASS_AND_RHS, GENERAL


class dynamic_step(sd_step):
    def __init__(self, model, number, name, previous, period, **kwds):
        super(dynamic_step, self).__init__(model, number, name, previous, period)
        for (key, val) in kwds.items():
            if key == "frames":
                key = "_frames"
            setattr(self, key, val)

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, period=1.0, frames=10, alpha=0.5, beta=0.0):

        period = getattr(self, "period", period)
        frames = getattr(self, "_frames", frames)

        time = np.array([0.0, self.start])
        dtime = period / float(frames)

        shape = (self.model.numdof, 2)
        un = np.zeros(shape)
        vn = np.zeros(shape)
        an = np.zeros(shape)

        # GET MASS AND RHS AT TIME 0
        Q = self.cload(time[0])
        dltyp, dlmag = self.dload(time[0])

        energy = ddlmag = mdload = None
        pnewdt = np.zeros(1)
        lflags = [DYNAMIC, SMALL_DISPLACEMENT, MASS_AND_RHS, GENERAL, 0]
        A = np.zeros((self.model.numdof, self.model.numdof))
        rhs = np.zeros(self.model.numdof)
        assemble_system(
            rhs,
            A,
            self.svtab,
            self.svars,
            energy,
            Q,
            self.dofs,
            un[:, 0],
            vn[:, 0],
            an[:, 0],
            time,
            dtime,
            0,
            1,
            1,
            dltyp,
            dlmag,
            self.predef,
            lflags,
            ddlmag,
            mdload,
            pnewdt,
            period,
            self.model.mesh.element_blocks,
            self.model.mesh.elemap,
            self.model.elements,
            self.model.eftab,
        )

        # LUMPED MASS MATRIX
        mass = np.array([sum(row) for row in A])
        an[:, 0] = rhs / mass

        for kframe in range(frames):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0] + dtime)
            X = self.dofvals(time[0] + dtime)
            dltyp, dload = self.dload(time[0] + dtime)

            rhs = np.zeros(self.model.numdof)
            A = np.zeros((self.model.numdof, self.model.numdof))
            assemble_system(
                rhs,
                A,
                self.svtab,
                self.svars,
                energy,
                Q,
                self.dofs,
                un[:, 0],
                vn[:, 0],
                an[:, 0],
                time,
                dtime,
                self.number,
                kframe + 1,
                1,
                dltyp,
                dlmag,
                self.predef,
                lflags,
                ddlmag,
                mdload,
                pnewdt,
                self.period,
                self.model.mesh.element_blocks,
                self.model.mesh.elemap,
                self.model.elements,
                self.model.eftab,
            )
            # LUMPED MASS MATRIX
            mass = np.array([sum(row) for row in A])

            # UPDATE ACCELERATION
            an[:, 1] = rhs / mass

            # UPDATE VELOCITY
            dv = ((1.0 - alpha) * an[:, 0] + alpha * an[:, 1]) * dtime
            vn[:, 1] = vn[:, 0] + dv

            # UPDATE DISPLACEMENT
            du = (
                vn[:, 0] + 0.5 * ((1.0 - beta) * an[:, 0] + beta * an[:, 1]) * dtime
            ) * dtime
            un[:, 1] = un[:, 0] + du
            un[self.doftags, 1] = X

            # ADVANCE KINEMATICS
            un[:, 0] = un[:, 1]
            vn[:, 0] = vn[:, 1]
            an[:, 0] = an[:, 1]

            time += dtime

            self.dofs[:] = un[:, 0]
            # self.vel[:] = vn[:,0]
            # self.accel[:] = an[:,0]

            self.advance(dtime, self.dofs)

        return
