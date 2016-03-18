from numpy import *
from copy import deepcopy

from .utilities import *
from .constants import *
from .data import *

class Step(object):
    def __init__(self, model, number, name, previous, period):
        self.model = model
        self.written = 0
        self.name = name
        if previous is None:
            self.time = 0.
        else:
            self.time = previous.frames[-1].value
        self.frames = []
        self.Frame(0.)
        self.period = period
        self.number = number
        self.previous = previous

        # DOF CONTAINERS
        self.dofs = zeros(self.model.numdof)
        self.dofx = {}

        # CONCENTRATED LOAD CONTAINERS
        self.cloadx = {}

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

    @property
    def doftags(self):
        return array(sorted(self.dofx), dtype=int)

    @property
    def dofvals(self):
        return array([self.dofx[key] for key in self.doftags])

    @property
    def cltags(self):
        return array(sorted(self.cloadx), dtype=int)

    @property
    def clvals(self):
        return array([self.cloadx[key] for key in self.cltags])

    def Frame(self, dtime, copy=1):
        frame = Frame(self.time, dtime)
        self.time += dtime
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = deepcopy(frame_n.field_outputs)
        frame.number = len(self.frames)
        self.frames.append(frame)
        return frame

    def copy_from(self, step):
        self.frames[0].field_outputs = deepcopy(step.frames[-1].field_outputs)
        self.dofs[:] = step.dofs
        self.dofx = dict(step.dofx)
        self.cloadx = dict(step.cloadx)
        self.dload[:] = step.dload
        self.dltyp[:] = step.dltyp
        self.predef[:] = step.predef
        self.svars[:] = step.svars

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
        self.assign_dof(DIRICHLET, nodes, ALL, 0.)
    FixDOF = FixNodes

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

        if dof == ALL:
            dofs = self.model.active_dof
        elif dof == PIN:
            dofs = [x for x in (X,Y,Z) if x in self.model.active_dof]
        elif not is_listlike(dof):
            dofs = [dof]
        else:
            dofs = dof

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

        for (i,inode) in enumerate(inodes):
            for j in dofs:
                I = self.model.dofmap(inode, j)
                if I is None:
                    raise UserInputError('INVALID DOF FOR NODE {0}'.format(inode))
                if I in self.cloadx and doftype == DIRICHLET:
                    msg = 'ATTEMPTING TO APPLY LOAD AND DISPLACEMENT '
                    msg += 'ON SAME DOF'
                    raise UserInputError(msg)
                elif I in self.dofx and doftype == NEUMANN:
                    msg = 'ATTEMPTING TO APPLY LOAD AND DISPLACEMENT '
                    msg += 'ON SAME DOF'
                    raise UserInputError(msg)
                if doftype == DIRICHLET:
                    self.dofx[I] = float(a[i])
                else:
                    self.cloadx[I] = float(a[i])

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def ConcentratedLoad(self, nodes, dof, amplitude=0.):
        self.assign_dof(NEUMANN, nodes, dof, amplitude)

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
