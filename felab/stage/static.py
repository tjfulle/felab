from numpy import *

from .stage import sd_stage
from ..constants import *
from ..utilities import *
from ..assembly import assemble_system, apply_boundary_conditions

class static_stage(sd_stage):
    procedure = STATIC
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

        elif solver == RIKS:
            self.riks_solve(**kwargs)

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
        K, rhs = assemble_system(
            self.model.mesh.element_blocks, self.model.mesh.elemap,
            self.model.elements, self.model.eftab, self.dofs, u, Qf, self.svtab,
            self.svars, dltyp, dload, self.predef, self.procedure, DIRECT,
            cflag=STIFF_AND_RHS, time=time)
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
        assemble_system(self.model.mesh.element_blocks, self.model.mesh.elemap,
                        self.model.elements, self.model.eftab, self.dofs, u, Qf,
                        self.svtab, self.svars, dltyp, dload, self.predef,
                        self.procedure, DIRECT, cflag=LP_OUTPUT)

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

        for iincrement in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for nit in range(maxit2):

                K, rhs = assemble_system(
                    self.model.mesh.element_blocks, self.model.mesh.elemap,
                    self.model.elements, self.model.eftab, self.dofs,
                    u, Q, self.svtab, self.svars, dltyp, dload, self.predef,
                    self.procedure, GENERAL, cflag=STIFF_AND_RHS, time=time,
                    istage=self.number, iincrement=iincrement+1, ninc=nit+1)

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

                if nit < maxit1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        logging.debug('CONVERGING TO LOSER TOLERANCE ON STAGE '
                                      '{0}, INCREMENT {1}'.format(self.number, iincrement+1))
                        break

                continue

            else:
                message  = 'FAILED TO CONVERGE ON STAGE '
                message += '{0}, INCREMENT {1}'.format(self.number, iincrement+1)
                logging.error(message)
                raise RuntimeError(message)

            logging.debug('STAGE {0}, INCREMENT {1}, '
                          'COMPLETE.'.format(self.number, iincrement+1))
            time += dtime
            self.dofs += u
            self.advance(dtime, self.dofs)

        return

    def riks_solve(self, period=1., increments=5, maxiters=20,
                   tolerance=1e-4, relax=1., tolerance1=1e-6):
        raise NotImplementedError("THIS FUNCTION IS STILL A NEWTON SOLVER")

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

        for iincrement in range(increments):

            # GET LOADS AND PRESCRIBED DISPLACEMENTS
            Q = self.cload(time[0]+dtime)
            X = self.dofvals(time[0]+dtime)
            dltyp, dload = self.dload(time[0]+dtime)

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for nit in range(maxit2):

                K, fext, fint = assemble_system(
                    self.model.mesh.element_blocks, self.model.mesh.elemap,
                    self.model.elements, self.model.eftab, self.dofs,
                    self.dofs, u, Q, self.svtab, self.svars, dltyp,
                    dload, self.predef, self.procedure, GENERAL,
                    cflag=STIFF_AND_RHS, time=time, istage=self.number,
                    iincrement=iincrement+1, ninc=nit+1, disp=1)
                rhs = fext - fint

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

                if nit < maxit1:
                    if err1 < tolerance1:
                        break
                else:
                    if err1 < tolerance:
                        break
                    elif err2 < 5e-2:
                        logging.debug('CONVERGING TO LOSER TOLERANCE ON STAGE '
                                      '{0}, INCREMENT {1}'.format(self.number, iincrement+1))
                        break

                continue

            else:
                message  = 'FAILED TO CONVERGE ON STAGE '
                message += '{0}, INCREMENT {1}'.format(self.number, iincrement+1)
                logging.error(message)
                raise RuntimeError(message)

            logging.debug('STAGE {0}, INCREMENT {1}, '
                          'COMPLETE.'.format(self.number, iincrement+1))
            time += dtime
            self.dofs += u
            self.advance(dtime, self.dofs)

        return
