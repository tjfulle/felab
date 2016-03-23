from numpy import *

from .constants import *
from .utilities import *
from .step import Step

class DynamicStep(Step):
    def __init__(self, model, number, name, previous, period, increments, nlgeom):
        super(DynamicStep, self).__init__(model, number, name, previous, period)
        self.nlgeom = nlgeom
        self.increments = increments

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
        self.dltyp[ielem].append(DLOAD)
        self.dload[ielem].append(asarray(components))

    def SurfaceLoad(self, surface, components):
        if not is_listlike(components):
            raise UserInputError('EXPECTED SURFACE LOAD VECTOR')
        if len(components) != self.model.dimensions:
            raise UserInputError('EXPECTED {0} SURFACE LOAD '
                                 'COMPONENTS'.format(len(self.model.active_dof)))
        surface = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surface:
            self.dltyp[iel].append(SLOAD)
            self.dload[iel].append([iedge]+[x for x in components])

    def SurfaceLoadN(self, surface, amplitude):
        surface = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surface:
            # DETERMINE THE NORMAL TO THE EDGE
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
    def run(self):

        period, increments = self.period, self.increments

        ti = self.frames[-1].start
        time = array([ti, ti])
        dtime = period / float(increments)

        shape = (self.model.numdof, 2)
        un = zeros(shape)
        vn = zeros(shape)
        an = zeros(shape)

        count = 0

        m = mass_matrix
        F = external_force_array_at_time_equal_0
        an[:,0] = F / m

        a, b = self.alpha, self.beta
        for n in range(increments):

            mass, rhs = self.model.assemble(self.dofs, u, Q,
                                            self.svtab, self.svars,
                                            self.dltyp, self.dload, self.predef,
                                            DYNAMIC, MASS_AND_RHS,
                                            istep=self.number, iframe=iframe+1)

            an[:,1] = rhs / mass
            vn[:,1] = vn[:,0] + ((1. - a) * an[:,0] + a * an[:,1]) * dtime
            un[:,1] = un[:,0] + (vn[:,0]+.5*((1.-b)*an[:,0]+b*an[:,1])*dtime)*dtime

            un[:,0] = un[:,1]
            vn[:,0] = vn[:,1]
            an[:,0] = an[:,1]

            time[1] += dtime

            self.dofs[:] = un[:,0]
            self.vel[:] = vn[:,0]
            self.accel[:] = an[:,0]
            self.advance(dtime, self.dofs)

        return
