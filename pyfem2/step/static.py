from numpy import *

from ..constants import *
from ..utilities import *
from ._step import SDStep

class StaticStep(SDStep):
    procedure = STATIC
    def __init__(self, model, number, name, previous, period=1., **kwds):

        super(StaticStep, self).__init__(model, number, name, previous, period)

        self.solver = kwds.pop('solver', DEFAULT)

        if self.solver == NEWTON:
            self.increments = kwds.pop('increments', 5)
            self.maxiters = kwds.pop('maxiters', 20)
            self.tolerance = kwds.pop('tolerance', 1e-4)
            self.relax = kwds.pop('relax', 1.)
            self.tolerance1 = kwds.pop('tolerance1', 1e-6)

        if self.solver not in (DEFAULT, NEWTON):
            raise TypeError('UNKNOWN STATIC SOLVER {0!r}'.format(self.solver))

        if kwds:
            keys = ', '.join([x for x in kwds.keys()])
            raise TypeError('UNKNOWN KEYWORDS FOR STATIC SOLVER {0!r}: '
                            '{1}'.format(self.solver, keys))

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self):

        if self.ran:
            raise RuntimeError('STEP ALREADY RUN')

        if self.solver == DEFAULT:
            self.direct_solve()

        elif self.solver == NEWTON:
            self.newton_solve()

        else:
            raise NotImplementedError

        self.ran = True

    def direct_solve(self):

        frame = self.Frame()
        time = array([0., self.start])

        # CONCENTRATED FORCES
        Qf = self.cload(self.period)
        dltyp, dload = self.dload(self.period)
        X = self.dofvals(self.period)

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        u = zeros_like(self.dofs)
        K, rhs = self.model.assemble(
            self.dofs, u, Qf, frame.field_outputs, dltyp, dload,
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
            self.dofs, u, Qf, frame.field_outputs, dltyp, dload,
            self.predef, self.procedure, DIRECT, cflag=LP_OUTPUT)

        self.dofs = u
        u, R, temp = self.model.format_dof(self.dofs)
        RF, M, Q = self.model.format_dof(react)
        frame.snapshot(self.period, U=u, R=R, T=temp, RF=RF, Q=Q, M=M)

    def newton_solve(self):

        period = self.period
        increments = self.increments
        maxiters = self.maxiters
        tolerance = self.tolerance
        relax = self.relax
        tolerance1 = self.tolerance1

        # TIME IS:
        # TIME[0]: VALUE OF STEP TIME AT BEGINNING OF INCREMENT
        # TIME[1]: VALUE OF TOTAL TIME AT BEGINNING OF INCREMENT
        time = array([0., self.start])
        dtime = period / float(increments)

        maxit2 = int(maxiters)
        maxit1 = max(int(maxit2/2.),1)

        for iframe in range(increments):

            frame = self.Frame()

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for nit in range(maxit2):

                K, rhs = self.model.assemble(
                    self.dofs, u, Q, frame.field_outputs, dltyp, dload,
                    self.predef, self.procedure, GENERAL, cflag=STIFF_AND_RHS,
                    time=time, istep=self.number, iframe=iframe+1, ninc=nit+1)

                # ENFORCE BOUNDARY CONDITIONS
                Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, X, self.dofs, u)

                # --- SOLVE FOR THE NODAL DISPLACEMENT
                w = linsolve(Kbc, Fbc)

                # --- UPDATE DISPLACEMENT INCREMENT
                u += relax * w

                # --- CHECK CONVERGENCE
                norm = lambda x_: sqrt(dot(x_, x_))
                err1 = sqrt(dot(w, w))
                dnom = sqrt(dot(u, u))
                if dnom > 1e-12:
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
            u, R, temp = self.model.format_dof(self.dofs)
            frame.snapshot(dtime, U=u, T=temp)

        return
