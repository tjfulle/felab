from numpy import *
from copy import deepcopy

from .utilities import *
from .constants import *
from .data import *

__all__ = ['StepRepository']

class StepRepository(object):
    def __init__(self, model):
        self.model = model
        self._keys = []
        self._values = []

    def __setitem__(self, key, value):
        self._keys.append(key)
        self._values.append(value)

    def __getitem__(self, key):
        return self._values[self._keys.index(key)]
        try:
            i = self._keys.index(key)
        except ValueError:
            raise KeyError(key)
        return self._values[i]

    def __len__(self):
        return len(self._keys)

    def __iter__(self):
        return iter(self._keys)

    def __contains__(self, key):
        return key in self._keys

    def keys(self):
        return self._keys

    def values(self):
        return self._values

    def items(self):
        for (i, key) in enumerate(self._keys):
            yield (key, self._values[i])

    def Step(self, name, type=GENERAL, period=1., increments=None, copy=1):
        if not self._values:
            step = Step(self.model, name, 0.)
        else:
            last = self._values[-1].frames[-1]
            if not last.converged:
                raise RuntimeError('PREVIOUS STEP HAS UNCONVERGED FRAMES')
            step = Step(self.model, name, last.value, type, period, increments)
            if copy:
                step.frames[0].field_outputs = deepcopy(last.field_outputs)
            step.frames[0].converged = True
        self[name] = step
        return self._values[-1]

    @property
    def last(self):
        return self._values[-1]

    @property
    def first(self):
        return self._values[0]

class Step(object):
    def __init__(self, model, name, start, type=None, period=None, increments=None):
        self.model = model
        self.written = 0
        self.name = name
        if type is None and abs(start) > 1e-12:
            raise ValueError('No step type given')
        self.time = start
        self.frames = []
        self.Frame(0.)

        self.type = type
        self.period = period
        self.increments = increments

        # DOF CONTAINERS
        self.dofs = zeros(self.model.numdof)
        self.doftags = []
        self.dofvals = []

        # CONCENTRATED LOAD CONTAINERS
        self.cltags = []
        self.clvals = []

        # CONTAINERS TO HOLD LOADS
        self.dload = emptywithlists(self.model.numele)
        self.dltyp = emptywithlists(self.model.numele)

        # PREDEFINED FIELDS
        self.predef = zeros((3, 1, self.model.numnod))

        # --- ALLOCATE STORAGE FOR SIMULATION DATA
        # STATE VARIABLE TABLE
        svtab = []
        nstatev = 0
        for el in self.model.elements:
            if not el.variables:
                svtab.append([])
                continue
            if not el.ndir:
                m = 1
            else:
                m = el.ndir + el.nshr
            m *= len(el.variables)
            if el.integration:
                m *= el.integration
            a = [nstatev, nstatev+m]
            svtab.append(slice(*a))
            nstatev += m
        self.svars = zeros((2, nstatev))
        self.svtab = svtab

    def __len__(self):
        return len(self.frames)

    def Frame(self, dtime, copy=1):
        frame = Frame(self.time, dtime)
        self.time += dtime
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = deepcopy(frame_n.field_outputs)
        frame.number = len(self.frames)
        self.frames.append(frame)
        return frame

    # ----------------------------------------------------------------------- #
    # --- BOUNDARY CONDITIONS ----------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def FixNodes(self, nodes):
        """Fix nodal degrees of freedom

        Parameters
        ----------
        nodes : int, list of int, or symbolic constant
            Nodes to fix

        Notes
        -----
        ``nodes`` can be a single external node label, a list of external node
        labels, or one of the region symbolic constants.

        All active displacement and rotation degrees of freedom are set to 0.

        """
        self.assign_dof(DIRICHLET, nodes, self.model.active_dof, 0.)
    FixDOF = FixNodes

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
        dof = [x for x in (X,Y,Z) if x in self.model.active_dof]
        self.assign_dof(DIRICHLET, nodes, dof, 0.)

    def PrescribedBC(self, nodes, dof, amplitude=0.):
        """Prescribe nodal degrees of freedom

        Parameters
        ----------
        nodes : int, list of int, or symbolic constant
            Nodes to fix
        dof : symbolic constant
            Degree of freedom to fix.  One of ``X,Y,Z,TX,TY,TZ,T``.
        amplitude : float or callable {0}
            The magnitude of the prescribed boundary condition

        Notes
        -----
        ``nodes`` can be a single external node label, a list of external node
        labels, or one of the region symbolic constants.

        ``amplitude`` can either be a float or a callable function. If a
        float, that value is assigned to all ``nodes``. If a callable
        function, the value assigned to each node is ``amplitude(x)``, where
        ``x`` is the node's coordinate position. The coordinate positions of
        all nodes are sent to the function as a n-dimensional column vector.

        Examples
        --------

        - Assign constant amplitude BC to the :math:`x` displacement of all
          nodes on left side of domain:

          .. code:: python

             self.PrescribedBC(ILO, X, 5.)

        - Assign variable amplitude BC to the :math:`x` displacement of all
          nodes on left side of domain. The variable amplitude function is
          :math:`\Delta_x=y^2`.

          .. code:: python

             fun = lambda x: x[:,1]**2
             self.PrescribedBC(ILO, X, fun)

        """
        self.assign_dof(DIRICHLET, nodes, dof, amplitude)
    PrescribedDOF = PrescribedBC

    def assign_dof(self, doftype, nodes, dof, amplitude):
        inodes = self.model.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, '__call__'):
            # AMPLITUDE IS A FUNCTION
            a = amplitude(self.model.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # CREATE A SINGLE AMPLITUDE FOR EACH NODE
            a = ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError('INCORRECT AMPLITUDE LENGTH')
            # AMPLITUDE IS A LIST OF AMPLITUDES
            a = asarray(amplitude)

        dofs = dof if is_listlike(dof) else [dof]
        for (i,inode) in enumerate(inodes):
            for j in dofs:
                I = self.model.dofmap(inode, j)
                if I is None:
                    raise UserInputError('INVALID DOF FOR NODE {0}'.format(inode))
                if I in self.cltags and doftype == DIRICHLET:
                    msg = 'ATTEMPTING TO APPLY LOAD AND DISPLACEMENT '
                    msg += 'ON SAME DOF'
                    raise UserInputError(msg)
                elif I in self.doftags and doftype == NEUMANN:
                    msg = 'ATTEMPTING TO APPLY LOAD AND DISPLACEMENT '
                    msg += 'ON SAME DOF'
                    raise UserInputError(msg)
                if doftype == DIRICHLET:
                    self.doftags.append(I)
                    self.dofvals.append(float(a[i]))
                else:
                    self.cltags.append(I)
                    self.clvals.append(float(a[i]))

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def ConcentratedLoad(self, nodes, dof, amplitude=0.):
        self.assign_dof(NEUMANN, nodes, dof, amplitude)

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
    # --- HEAT TRANSFER LOADINGS -------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def SurfaceFlux(self, surface, qn):
        surf = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.dltyp[iel].append(SFLUX)
            self.dload[iel].append([iedge, qn])

    def SurfaceConvection(self, surface, Too, h):
        if self.model.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        surf = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.dltyp[iel].append(SFILM)
            self.dload[iel].append([iedge, Too, h])

    def HeatGeneration(self, region, amplitude):
        if region == ALL:
            xelems = sorted(self.model.mesh.elemap.keys())
        else:
            xelems = region
        inodes = []
        for eb in self.model.mesh.eleblx:
            for (i, xel) in enumerate(eb.labels):
                if xel in xelems:
                    inodes.extend(eb.elecon[i])
        inodes = unique(inodes)
        if hasattr(amplitude, '__call__'):
            # AMPLITUDE IS A FUNCTION
            x = self.model.mesh.coord[inodes]
            a = amplitude(x)
        elif not is_listlike(amplitude):
            a = amplitude * ones(len(inodes))
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError('HEAT GENERATION AMPLITUDE MUST HAVE '
                                     'LENGTH {0}'.format(len(inodes)))
            a = asarray(amplitude)
        nodmap = dict(zip(inodes, range(inodes.shape[0])))
        for xelem in xelems:
            ielem = self.model.mesh.elemap[xelem]
            ix = [nodmap[n] for n in self.model.elements[ielem].inodes]
            self.dltyp[ielem].append(HSRC)
            self.dload[ielem].append(a[ix])

    def Temperature(self, nodes, amplitude):
        inodes = self.model.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, '__call__'):
            # AMPLITUDE IS A FUNCTION
            a = amplitude(self.model.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # CREATE A SINGLE AMPLITUDE FOR EACH NODE
            a = ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError('INCORRECT AMPLITUDE LENGTH')
            # AMPLITUDE IS A LIST OF AMPLITUDES
            a = asarray(amplitude)
        self.final_temp = a

    # ----------------------------------------------------------------------- #
    # --- RUN AND SOLVING --------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def run(self):
        if self.type == LINEAR_PERTURBATION:
            return self._linear_perturbation_step()
        elif self.type == GENERAL:
            return self._general_static_step()
        raise NotImplementedError

    def _linear_perturbation_step(self):

        # CONCENTRATED FORCES
        Q = zeros_like(self.dofs)
        Q[self.cltags] = self.clvals

        # ASSEMBLE THE GLOBAL STIFFNESS AND FORCE
        du = zeros_like(self.dofs)
        K, rhs = self.model.assemble(self.dofs, du, Q, self.svtab, self.svars,
                                     self.dltyp, self.dload, self.predef,
                                     step_type=LINEAR_PERTURBATION)

        # APPLY BOUNDARY CONDITIONS
        Kbc, Fbc = self.model.apply_bc(K, rhs, self.doftags, self.dofvals)

        # SOLVE FOR UNKNOWN DOFS
        u = linsolve(Kbc, Fbc)

        # SANITY CHECK
        if not allclose(u[self.doftags], self.dofvals):
            logging.warn('INCORRECT SOLUTION TO DOFS')

        # TOTAL FORCE, INCLUDING REACTION, AND REACTION
        R = dot(K, u) - rhs
        self.model.assemble(self.dofs, u, Q, self.svtab, self.svars,
                            self.dltyp, self.dload, self.predef,
                            step_type=LINEAR_PERTURBATION, cflag=LP_OUTPUT)
        self.dofs += u
        self.advance(self.period, R=R)

    def advance(self, dtime, **kwds):
        frame_n = self.frames[-1]
        if not frame_n.converged:
            raise RuntimeError('ATTEMPTING TO UPDATE AN UNCONVERGED FRAME')

        # ADVANCE STATE VARIABLES
        self.svars[0] = self.svars[1]

        # CREATE FRAME TO HOLD RESULTS
        frame = self.Frame(dtime)

        # STORE DEGREES OF FREEDOM
        u = self.dofs.reshape(self.model.mesh.coord.shape)
        frame.field_outputs['U'].add_data(u)

        # STORE KEYWORDS
        for (kwd, val) in kwds.items():
            frame.field_outputs[kwd].add_data(val)

        for (ieb, eb) in enumerate(self.model.mesh.eleblx):
            if not eb.eletyp.variables:
                continue

            # PASS VALUES FROM SVARS TO THE FRAME OUTPUT
            ntens = eb.eletyp.ndir + eb.eletyp.nshr
            m = 1 if not eb.eletyp.integration else eb.eletyp.integration
            n = len(eb.eletyp.variables)
            for (e, xel) in enumerate(eb.labels):
                iel = self.model.mesh.elemap[xel]
                el = self.model.elements[iel]
                ue = u[el.inodes]
                svars = self.svars[0,self.svtab[iel]].reshape(m,n,ntens)
                for (j, name) in enumerate(el.variables):
                    frame.field_outputs[eb.name,name].add_data(svars[:,j], e)

        frame.converged = True

class Frame(object):
    def __init__(self, start, dtime):
        self.start = start
        self.increment = dtime
        self.value = start + dtime
        self.field_outputs = FieldOutputs()
        self.converged = False

    def adjust_dt(self, dtime):
        self.increment = dtime
        self.value = self.start + dtime

    def SymmetricTensorField(self, name, position, labels, ndir, nshr,
                             eleblk=None, ngauss=None,  elements=None, data=None):
        field = SymmetricTensorField(name, position, labels, ndir, nshr,
                                     eleblk=eleblk, ngauss=ngauss,
                                     elements=elements, data=data)

        if field.eleblk is not None:
            name = (field.eleblk, name)
        self.field_outputs[name] = field

        return field

    def VectorField(self, name, position, labels, ncomp, eleblk=None,
                    ngauss=None, elements=None, data=None):

        field = VectorField(name, position, labels, ncomp, eleblk=eleblk,
                            ngauss=ngauss, elements=elements, data=data)

        if field.eleblk is not None:
            name = (field.eleblk, name)
        self.field_outputs[name] = field

        return field

    def ScalarField(self, name, position, labels, eleblk=None, ngauss=None,
                    elements=None, data=None):

        field = ScalarField(name, position, labels, eleblk=eleblk,
                            ngauss=ngauss, elements=elements, data=data)

        if field.eleblk is not None:
            name = (field.eleblk, name)
        self.field_outputs[name] = field

        return field

    def add_data(self, **kwds):
        for (key, value) in kwds.items():
            d = {}
            if isinstance(value, tuple):
                if len(value) == 2:
                    d['ix'] = value[1]
                    value = value[0]
                else:
                    raise ValueError('Unknown add_data option for {0}'.format(key))
            self.field_outputs[key].add_data(value, **d)
