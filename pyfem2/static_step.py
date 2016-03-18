from numpy import *

from .constants import *
from .utilities import *
from .step import Step

class StaticStep(Step):
    def __init__(self, model, number, name, previous, period, increments,
                 maxiters, nlgeom, solver):
        super(StaticStep, self).__init__(model, number, name, previous, period)
        self.nlgeom = nlgeom
        self.solver = solver
        self.increments = increments
        self.maxiters = maxiters

    def PinNodes(self, nodes):
        """Pin nodal degrees of freedom

        Parameters
        ----------
        nodes : int, list of int, or symbolic constant
            Nodes to fix

        Notes
        -----
        ``nodes`` can be a single external node label, a list of external node
        labels, or one of the region symbolic constants.

        All active displacement degrees of freedom are set to 0.

        """
        self.assign_dof(DIRICHLET, nodes, PIN, 0.)

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def GravityLoad(self, region, components):
        orphans = self.model.orphaned_elements
        if region == ALL:
            ielems = range(self.model.numele)
        else:
            ielems = [self.model.mesh.elemap[el] for el in region]
        if not is_listlike(components):
            raise UserInputError('EXPECTED GRAVITY LOAD VECTOR')
        if len(components) != self.model.dimensions:
            raise UserInputError('EXPECTED {0} GRAVITY LOAD '
                                 'COMPONENTS'.format(len(self.model.active_dof)))
        a = asarray(components)
        for iel in ielems:
            if iel in orphans:
                raise UserInputError('ELEMENT BLOCKS MUST BE ASSIGNED '
                                     'BEFORE GRAVITY LOADS')
            el = self.model.elements[iel]
            rho = el.material.density
            if rho is None:
                raise UserInputError('ELEMENT MATERIAL DENSITY MUST BE ASSIGNED '
                                     'BEFORE GRAVITY LOADS')
            self.dltyp[iel].append(DLOAD)
            self.dload[iel].append(rho*a)

    def DistributedLoad(self, region, components):
        if not is_listlike(components):
            raise UserInputError('EXPECTED DISTRIBUTED LOAD VECTOR')
        if len(components) != self.model.dimensions:
            raise UserInputError('EXPECTED {0} DISTRIBUTED LOAD '
                                 'COMPONENTS'.format(len(self.model.active_dof)))
        if region == ALL:
            ielems = range(self.model.numele)
        elif not is_listlike(region):
            ielems = [self.model.mesh.elemap[region]]
        else:
            ielems = [self.model.mesh.elemap[el] for el in region]
        if any(in1d(ielems, self.model.orphaned_elements)):
            raise UserInputError('ELEMENT PROPERTIES MUST BE ASSIGNED '
                                 'BEFORE DISTRIBUTEDLOAD')
        self.dltyp[ielem].append(DLOAD)
        self.dload[ielem].append(asarray(components))

    def SurfaceLoad(self, surface, components):
        if not is_listlike(components):
            raise UserInputError('EXPECTED SURFACE LOAD VECTOR')
        if len(components) != self.model.dimensions:
            raise UserInputError('EXPECTED {0} SURFACE LOAD '
                                 'COMPONENTS'.format(len(self.model.active_dof)))
        surface = self.model.mesh.find_surface(surface)
        orphans = self.model.orphaned_elements
        for (iel, iedge) in surface:
            if iel in orphans:
                raise UserInputError('ELEMENT PROPERTIES MUST BE ASSIGNED '
                                     'BEFORE SURFACELOAD')
            self.dltyp[iel].append(SLOAD)
            self.dload[iel].append([iedge]+[x for x in components])

    def SurfaceLoadN(self, surface, amplitude):
        surface = self.model.mesh.find_surface(surface)
        orphans = self.model.orphaned_elements
        for (iel, iedge) in surface:
            # DETERMINE THE NORMAL TO THE EDGE
            if iel in orphans:
                raise UserInputError('ELEMENT PROPERTIES MUST BE ASSIGNED '
                                     'BEFORE SURFACELOADN')
            el = self.model.elements[iel]
            edgenod = el.edges[iedge]
            xb = el.xc[edgenod]
            if self.model.dimensions == 2:
                n = normal2d(xb)
            else:
                raise NotImplementedError('3D SURFACE NORMAL')
            self.dltyp[iel].append(SLOAD)
            self.dload[iel].append([iedge]+[x for x in amplitude*n])

    def Pressure(self, surface, amplitude):
        self.SurfaceLoadN(surface, -amplitude)

    # ----------------------------------------------------------------------- #
    # --- RUN --------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self, **kwargs):

        if self.solver in (NEWTON, ):
            if self.increments is None:
                self.increments = 5

        if self.increments is None:
            return self.direct_solve()

        if self.solver in (NEWTON, None):
            return self.newton_solve(**kwargs)

        raise NotImplementedError

    def direct_solve(self):

        # CONCENTRATED FORCES
        Qf = zeros_like(self.dofs)
        Qf[self.cltags] = self.clvals

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        du = zeros_like(self.dofs)
        K, rhs = self.model.assemble(self.dofs, du, Qf, self.svtab, self.svars,
                                     self.dltyp, self.dload, self.predef,
                                     STATIC, DIRECT)

        # APPLY BOUNDARY CONDITIONS
        Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, self.dofvals)

        # SOLVE FOR UNKNOWN DOFS
        u = linsolve(Kbc, Fbc)

        # SANITY CHECK
        if not allclose(u[self.doftags], self.dofvals):
            logging.warn('INCORRECT SOLUTION TO DOFS')

        # TOTAL FORCE, INCLUDING REACTION, AND REACTION
        R = dot(K, u) - rhs
        self.model.assemble(self.dofs, u, Qf, self.svtab, self.svars,
                            self.dltyp, self.dload, self.predef,
                            STATIC, DIRECT, cflag=LP_OUTPUT)
        self.dofs += u
        self.advance(self.period, R=R)

    def newton_solve(self, tolerance=1e-4, relax=1., tolerance1=1e-6):

        period, increments, maxiters = self.period, self.increments, self.maxiters

        ti = self.frames[-1].start
        time = array([ti, ti])
        dtime = period / float(increments)

        maxit2 = int(maxiters)
        maxit1 = max(int(maxit2/2.),1)

        # CONCENTRATED FORCES AT END OF LAST STEP AND END OF THIS STEP
        Q0, Qf = zeros_like(self.dofs), zeros_like(self.dofs)
        Q0[self.previous.cltags] = self.previous.clvals
        Qf[self.cltags] = self.clvals

        # DISPLACEMENT BOUNDARY CONDITIONS AT END OF LAST STEP AND END OF THIS
        # STEP
        Xf = self.dofvals
        X0 = array([self.previous.dofx.get(I, 0) for I in self.doftags])

        for iframe in range(increments):

            load_fac = float(iframe+1) / float(increments)

            # INTERPOLATED CONCENTRATED LOADS AND DISPLACEMENTS
            Q = (1. - load_fac) * Q0 + load_fac * Qf
            X = (1. - load_fac) * X0 + load_fac * Xf

            # NEWTON-RAPHSON LOOP
            err1 = 1.
            u = zeros(self.model.numdof)
            for nit in range(maxit2):

                K, rhs = self.model.assemble(self.dofs, u, Q,
                                             self.svtab, self.svars,
                                             self.dltyp, self.dload, self.predef,
                                             STATIC, GENERAL,
                                             istep=self.number, iframe=iframe+1,
                                             ninc=nit+1, load_fac=load_fac,)
                # ENFORCE DISPLACEMENT BOUNDARY CONDITIONS
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
            time[1] += dtime
            self.dofs += u
            self.advance(dtime)

        return
