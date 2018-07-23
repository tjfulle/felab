from numpy import *

from .stage import sd_stage
from ..x.constants import *
from ..x.utilities import *
from ..assembly import assemble_system, apply_boundary_conditions

class static_stage(sd_stage):
    def __init__(self, model, number, name, previous, period=1., **kwds):
        super(static_stage, self).__init__(model, number, name, previous, period)
        for (key, val) in kwds.items():
            if key == 'increments':
                key = '_increments'
            setattr(self, key, val)

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, **kwargs):

        if self.ran:
            raise RuntimeError('STAGE ALREADY RUN')

        solver = kwargs.pop('solver', getattr(self, 'solver', None))
        increments = kwargs.get('increments', getattr(self, '_increments', None))

        if solver is None and increments:
            solver = NEWTON

        if solver is None:
            self.direct_solve()

        elif solver == NEWTON:
            self.newton_solve(**kwargs)

        else:
            raise NotImplementedError

        self.ran = True

    def direct_solve(self):

        time = array([0., self.start])

        # CONCENTRATED FORCES
        Qf = self.cload(self.period)
        dltyp, dlmag = self.dload(self.period)
        X = self.dofvals(self.period)

        # define arguments needed for assembly
        u = zeros(self.model.numdof)
        energy = v = a = ddlmag = mdload = pnewdt = None
        time = array([0., self.period])
        dtime = 1.
        lflags = [STATIC_DIRECT, SMALL_DISPLACEMENT, STIFF_AND_RHS, GENERAL, 0]

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        rhs = zeros(self.model.numdof)
        K = zeros((self.model.numdof, self.model.numdof))
        assemble_system(rhs, K, self.svtab, self.svars, energy, Qf,
                        self.dofs, u, v, a, time, dtime, self.number, 1, 0,
                        dltyp, dlmag, self.predef, lflags,
                        ddlmag, mdload, pnewdt, self.period,
                        self.model.mesh.element_blocks, self.model.mesh.elemap,
                        self.model.elements, self.model.eftab)
        self._K = K

        # ENFORCE BOUNDARY CONDITIONS
        Kbc, Fbc = apply_boundary_conditions(K, rhs, self.doftags, X)

        # SOLVE FOR UNKNOWN DOFS
        u[:] = linsolve(Kbc, Fbc)

        # SANITY CHECK
        if not allclose(u[self.doftags], X):
            logging.warn('INCORRECT SOLUTION TO DOFS')

        # TOTAL FORCE, INCLUDING REACTION, AND REACTION
        react = dot(K, u) - rhs

        # ASSEMBLE AGAIN - ONLY TO UPDATE STRESS IN ELEMENTS TO COMPUTED
        # DISPLACEMENT
        rhs2 = zeros_like(rhs)
        K2 = zeros_like(K)
        assemble_system(rhs2, K2, self.svtab, self.svars, energy, Qf,
                        self.dofs, u, v, a, time, dtime, self.number, 0, 0,
                        dltyp, dlmag, self.predef, lflags,
                        ddlmag, mdload, pnewdt, self.period,
                        self.model.mesh.element_blocks, self.model.mesh.elemap,
                        self.model.elements, self.model.eftab)

        self.dofs = u
        self.advance(self.period, self.dofs, react=react)

    def newton_solve(self, period=1., increments=5, maxiters=20,
                     tolerance=1e-4, relax=1., tolerance1=1e-6):

        period = getattr(self, 'period', period)
        increments = getattr(self, '_increments', increments)
        maxiters = getattr(self, 'maxiters', maxiters)

        # TIME IS:
        # TIME[0]: VALUE OF STAGE TIME AT BEGINNING OF INCREMENT
        # TIME[1]: VALUE OF TOTAL TIME AT BEGINNING OF INCREMENT
        time = array([0., self.start])
        dtime = period / float(increments)

        maxit2 = int(maxiters)
        maxit1 = max(int(maxit2/2.),1)

        energy = v = a = ddlmag = mdload = pnewdt = None
        lflags = [STATIC_ITERATIVE, SMALL_DISPLACEMENT, STIFF_AND_RHS, GENERAL, 0]
        for kinc in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dlmag = self.dload(time[0]+dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for kiter in range(maxit2):

                rhs = zeros(self.model.numdof)
                K = zeros((self.model.numdof, self.model.numdof))
                assemble_system(rhs, K, self.svtab, self.svars, energy, Q,
                                self.dofs, u, v, a, time, dtime, self.number,
                                kinc+1, kiter+1, dltyp, dlmag, self.predef, lflags,
                                ddlmag, mdload, pnewdt, self.period,
                                self.model.mesh.element_blocks,
                                self.model.mesh.elemap, self.model.elements,
                                self.model.eftab)

                # ENFORCE BOUNDARY CONDITIONS
                Kbc, Fbc = apply_boundary_conditions(K, rhs, self.doftags, X, self.dofs, u)

                # --- SOLVE FOR THE NODAL DISPLACEMENT
                w = linsolve(Kbc, Fbc)

                # --- UPDATE DISPLACEMENT INCREMENT
                u += relax * w

                # --- CHECK CONVERGENCE
                err1 = sqrt(dot(w, w))
                dnom = sqrt(dot(u, u))
                if dnom > 1e-8:
                    err1 /= dnom
                err2 = sqrt(dot(rhs, rhs)) / float(self.model.numdof)

                if kiter < maxit1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        logging.debug('CONVERGING TO LOSER TOLERANCE ON STAGE '
                                      '{0}, INCREMENT {1}'.format(self.number,
                                                                  kinc+1))
                        break

                continue

            else:
                message  = 'FAILED TO CONVERGE ON STAGE '
                message += '{0}, INCREMENT {1}'.format(self.number, kinc+1)
                logging.error(message)
                raise RuntimeError(message)

            logging.debug('STAGE {0}, INCREMENT {1}, '
                          'COMPLETE.'.format(self.number, kinc+1))
            time += dtime
            self.dofs += u
            self.advance(dtime, self.dofs)

        return
