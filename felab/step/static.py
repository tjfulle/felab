import numpy as np

import felab.util.tty as tty
from felab.constants import (
    NEWTON,
    GENERAL,
    STIFF_AND_RHS,
    SMALL_DISPLACEMENT,
    STATIC_DIRECT,
    STATIC_ITERATIVE,
)
from felab.util.numeric import linsolve
from felab.step.step import StressDisplacmentStep
from felab.assembly import assemble_system, apply_boundary_conditions


class StaticStep(StressDisplacmentStep):
    def __init__(self, model, number, name, previous, period=1.0, **kwds):
        super(StaticStep, self).__init__(model, number, name, previous, period)
        for (key, val) in kwds.items():
            if key == "frames":
                key = "_frames"
            setattr(self, key, val)

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, **kwargs):

        if self.ran:
            raise RuntimeError("STEP ALREADY RUN")

        solver = kwargs.pop("solver", getattr(self, "solver", None))
        frames = kwargs.get("frames", getattr(self, "_frames", None))

        if solver is None and frames:
            solver = NEWTON

        if solver is None:
            self.direct_solve()

        elif solver == NEWTON:
            self.newton_solve(**kwargs)

        else:
            raise NotImplementedError

        self.ran = True

    def direct_solve(self):

        time = np.array([0.0, self.start])

        # CONCENTRATED FORCES
        Qf = self.cload(self.period)
        dltyp, dlmag = self.dload(self.period)
        X = self.dofvals(self.period)

        # define arguments needed for assembly
        u = np.zeros(self.model.numdof)
        energy = v = a = ddlmag = mdload = pnewdt = None
        time = np.array([0.0, self.period])
        dtime = 1.0
        lflags = [STATIC_DIRECT, SMALL_DISPLACEMENT, STIFF_AND_RHS, GENERAL, 0]

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        rhs = np.zeros(self.model.numdof)
        K = np.zeros((self.model.numdof, self.model.numdof))
        assemble_system(
            rhs,
            K,
            self.svtab,
            self.svars,
            energy,
            Qf,
            self.dofs,
            u,
            v,
            a,
            time,
            dtime,
            self.number,
            1,
            0,
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
        self._K = K

        # ENFORCE BOUNDARY CONDITIONS
        Kbc, Fbc = apply_boundary_conditions(K, rhs, self.doftags, X)

        # SOLVE FOR UNKNOWN DOFS
        u[:] = linsolve(Kbc, Fbc)

        # SANITY CHECK
        if not np.allclose(u[self.doftags], X):
            tty.warn("INCORRECT SOLUTION TO DOFS")

        # TOTAL FORCE, INCLUDING REACTION, AND REACTION
        react = np.dot(K, u) - rhs

        # ASSEMBLE AGAIN - ONLY TO UPDATE STRESS IN ELEMENTS TO COMPUTED
        # DISPLACEMENT
        rhs2 = np.zeros_like(rhs)
        K2 = np.zeros_like(K)
        assemble_system(
            rhs2,
            K2,
            self.svtab,
            self.svars,
            energy,
            Qf,
            self.dofs,
            u,
            v,
            a,
            time,
            dtime,
            self.number,
            0,
            0,
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

        self.dofs = u
        self.advance(self.period, self.dofs, react=react)

    def newton_solve(
        self,
        period=1.0,
        frames=5,
        maxiters=20,
        tolerance=1e-4,
        relax=1.0,
        tolerance1=1e-6,
    ):

        period = getattr(self, "period", period)
        frames = getattr(self, "_frames", frames)
        maxiters = getattr(self, "maxiters", maxiters)

        # TIME IS:
        # TIME[0]: VALUE OF STEP TIME AT BEGINNING OF INCREMENT
        # TIME[1]: VALUE OF TOTAL TIME AT BEGINNING OF INCREMENT
        time = np.array([0.0, self.start])
        dtime = period / float(frames)

        maxit2 = int(maxiters)
        maxit1 = max(int(maxit2 / 2.0), 1)

        energy = v = a = ddlmag = mdload = pnewdt = None
        lflags = [STATIC_ITERATIVE, SMALL_DISPLACEMENT, STIFF_AND_RHS, GENERAL, 0]
        for kframe in range(frames):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0] + dtime)
            X = self.dofvals(time[0] + dtime)
            dltyp, dlmag = self.dload(time[0] + dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.0
            u = np.zeros(self.model.numdof)
            for kiter in range(maxit2):

                rhs = np.zeros(self.model.numdof)
                K = np.zeros((self.model.numdof, self.model.numdof))
                assemble_system(
                    rhs,
                    K,
                    self.svtab,
                    self.svars,
                    energy,
                    Q,
                    self.dofs,
                    u,
                    v,
                    a,
                    time,
                    dtime,
                    self.number,
                    kframe + 1,
                    kiter + 1,
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

                # ENFORCE BOUNDARY CONDITIONS
                Kbc, Fbc = apply_boundary_conditions(
                    K, rhs, self.doftags, X, self.dofs, u
                )

                # --- SOLVE FOR THE NODAL DISPLACEMENT
                w = linsolve(Kbc, Fbc)

                # --- UPDATE DISPLACEMENT INCREMENT
                u += relax * w

                # --- CHECK CONVERGENCE
                err1 = np.sqrt(np.dot(w, w))
                dnom = np.sqrt(np.dot(u, u))
                if dnom > 1e-8:
                    err1 /= dnom
                err2 = np.sqrt(np.dot(rhs, rhs)) / float(self.model.numdof)

                if kiter < maxit1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        tty.debug(
                            "CONVERGING TO LOSER TOLERANCE ON STEP "
                            "{0}, FRAME {1}".format(self.number, kframe + 1)
                        )
                        break

                continue

            else:
                message = "FAILED TO CONVERGE ON STEP "
                message += "{0}, FRAME {1}".format(self.number, kframe + 1)
                tty.error(message)
                raise RuntimeError(message)

            tty.debug(
                "STEP {0}, FRAME {1}, " "COMPLETE.".format(self.number, kframe + 1)
            )
            time += dtime
            self.dofs += u
            self.advance(dtime, self.dofs)

        return
