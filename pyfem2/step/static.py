from numpy import *

from ..constants import *
from ..utilities import *
from ._step import SDStep

class StaticStep(SDStep):
    procedure = STATIC
    def __init__(self, model, number, name, previous, period=1., **kwds):
        super(StaticStep, self).__init__(model, number, name, previous, period)
        for (key, val) in kwds.items():
            setattr(self, key, val)

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, **kwargs):

        if self.ran:
            raise RuntimeError('STEP ALREADY RUN')

        solver = kwargs.pop('solver', getattr(self, 'solver', None))
        increments = kwargs.get('increments', getattr(self, 'increments', None))

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
        dltyp, dload = self.dload(self.period)
        X = self.dofvals(self.period)

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        u = zeros_like(self.dofs)
        K, rhs = self.model.assemble(
            self.dofs, u, Qf, self.svtab, self.svars, dltyp, dload,
            self.predef, self.procedure, DIRECT, cflag=STIFF_AND_RHS, time=time)

        # ENFORCE BOUNDARY CONDITIONS
        Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, X)

        # SOLVE FOR UNKNOWN DOFS
        u[:] = linsolve(Kbc, Fbc)

        # SANITY CHECK
        if not allclose(u[self.doftags], X):
            logging.warn('INCORRECT SOLUTION TO DOFS')

        # TOTAL FORCE, INCLUDING REACTION, AND REACTION
        react = dot(K, u) - rhs

        # ASSEMBLE AGAIN - ONLY TO UPDATE STRESS IN ELEMENTS TO COMPUTED
        # DISPLACEMENT
        self.model.assemble(
            self.dofs, u, Qf, self.svtab, self.svars, dltyp, dload,
            self.predef, self.procedure, DIRECT, cflag=LP_OUTPUT)

        self.dofs = u
        self.advance(self.period, self.dofs, react=react)

    def newton_solve(self, period=1., increments=5, maxiters=20,
                     tolerance=1e-4, relax=1., tolerance1=1e-6):

        period = getattr(self, 'period', period)
        increments = getattr(self, 'increments', increments)
        maxiters = getattr(self, 'maxiters', maxiters)

        # TIME IS:
        # TIME[0]: VALUE OF STEP TIME AT BEGINNING OF INCREMENT
        # TIME[1]: VALUE OF TOTAL TIME AT BEGINNING OF INCREMENT
        time = array([0., self.start])
        dtime = period / float(increments)

        maxit2 = int(maxiters)
        maxit1 = max(int(maxit2/2.),1)

        for iframe in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for nit in range(maxit2):

                K, rhs = self.model.assemble(
                    self.dofs, u, Q, self.svtab, self.svars, dltyp, dload,
                    self.predef, self.procedure, GENERAL, cflag=STIFF_AND_RHS,
                    time=time, istep=self.number, iframe=iframe+1, ninc=nit+1)

                # ENFORCE BOUNDARY CONDITIONS
                Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, X, self.dofs, u)

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

                if nit < maxit1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        logging.debug('CONVERGING TO LOSER TOLERANCE ON STEP '
                                      '{0}, FRAME {1}'.format(self.number, iframe+1))
                        break

                continue

            else:
                message  = 'FAILED TO CONVERGE ON STEP '
                message += '{0}, FRAME {1}'.format(self.number, iframe+1)
                logging.error(message)
                raise RuntimeError(message)

            logging.debug('STEP {0}, FRAME {1}, '
                          'COMPLETE.'.format(self.number, iframe+1))
            time += dtime
            self.dofs += u
            self.advance(dtime, self.dofs)

        return

    def newton_solve1(self, period=1., increments=5, maxiters=20,
                      tolerance=1e-4, relax=1., tolerance1=1e-6):

        period = getattr(self, 'period', period)
        increments = getattr(self, 'increments', increments)
        maxiters = getattr(self, 'maxiters', maxiters)

        # TIME IS:
        # TIME[0]: VALUE OF STEP TIME AT BEGINNING OF INCREMENT
        # TIME[1]: VALUE OF TOTAL TIME AT BEGINNING OF INCREMENT
        time = array([0., self.start])
        dtime = period / float(increments)

        maxit2 = int(maxiters)
        maxit1 = max(int(maxit2/2.),1)

        for iframe in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for nit in range(maxit2):

                K, fext, fint = self.model.assemble(
                    self.dofs, u, Q, self.svtab, self.svars, dltyp, dload,
                    self.predef, self.procedure, GENERAL, cflag=STIFF_AND_RHS,
                    time=time, istep=self.number, iframe=iframe+1, ninc=nit+1,
                    disp=1)
                rhs = fext - fint

                # ENFORCE BOUNDARY CONDITIONS
                Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, X, self.dofs, u)

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

                if nit < maxit1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        logging.debug('CONVERGING TO LOSER TOLERANCE ON STEP '
                                      '{0}, FRAME {1}'.format(self.number, iframe+1))
                        break

                continue

            else:
                message  = 'FAILED TO CONVERGE ON STEP '
                message += '{0}, FRAME {1}'.format(self.number, iframe+1)
                logging.error(message)
                raise RuntimeError(message)

            logging.debug('STEP {0}, FRAME {1}, '
                          'COMPLETE.'.format(self.number, iframe+1))
            time += dtime
            self.dofs += u
            self.advance(dtime, self.dofs)

        return
