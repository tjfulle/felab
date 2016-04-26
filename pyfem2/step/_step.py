from numpy import *
from copy import deepcopy

from ..utilities import *
from ..constants import *
from .data_wharehouse import *

class Step(object):
    def __init__(self, model, number, name, previous, period):
        self.model = model
        self.written = 0
        self.ran = False
        self.name = name
        if previous is None:
            self.start = 0.
        else:
            self.start = previous.frames[-1].value
        self.frames = []
        self.Frame()
        self.period = period
        self.number = number
        self.previous = previous

        self.dofs = zeros(self.model.numdof)

        # DOFX[I] IS THE PRESCRIBED DOF FOR DOF I
        self.dofx = {}

        # CLOADX[I] IS THE PRESCRIBED CONCENTRATED LOAD FOR DOF I
        self.cloadx = {}

        # CONTAINERS TO HOLD DISTRIBUTED LOADS
        self.dloadx = {}
        self.sloadx = {}

        # CONTAINERS TO HOLD HEAT TRANSFER LOADS
        self.sfluxx = {}
        self.sfilmx = {}
        self.hsrcx = {}

        # PREDEFINED FIELDS
        self.predef = zeros((3, 1, self.model.numnod))

    def __len__(self):
        return len(self.frames)

    @property
    def doftags(self):
        return array(sorted(self.dofx), dtype=int)

    def dofvals(self, step_time):

        ix = self.doftags

        # DOFS AT END OF LAST STEP
        X0 = array([self.previous.dofx.get(I, 0) for I in ix])

        # DOFS AT END OF THIS STEP
        Xf = array([self.dofx[I] for I in ix])

        # INTERPOLATE CONCENTRATED LOAD TO CURRENT TIME
        fac = min(1., step_time / self.period)
        return (1. - fac) * X0 + fac * Xf

    @property
    def cltags(self):
        return array(sorted(self.cloadx), dtype=int)

    def cload(self, step_time):
        # CONCENTRATED LOAD AT END OF LAST STEP
        ix = self.previous.cltags
        Q0 = zeros_like(self.dofs)
        Q0[ix] = [self.previous.cloadx[key] for key in ix]

        # CONCENTRATED LOAD AT END OF THIS STEP
        ix = self.cltags
        Qf = zeros_like(self.dofs)
        Qf[ix] = [self.cloadx[key] for key in ix]

        # INTERPOLATE CONCENTRATED LOAD TO CURRENT TIME
        fac = min(1., step_time / self.period)
        return (1. - fac) * Q0 + fac * Qf

    def dload(self, step_time):

        # INTERPOLATES ALL DISTRIBUTED LOADS (BODY AND SURFACE) TO STEP_TIME

        # CONTAINER FOR ALL DLOADS
        dltyp = emptywithlists(self.model.numele)
        dload = emptywithlists(self.model.numele)

        # INTERPOLATION FACTOR
        fac = min(1., step_time / self.period)

        # INTERPOLATE SURFACE LOADS
        for (key, Ff) in self.sloadx.items():
            iel, iedge = key
            F0 = self.previous.sloadx.get(key, zeros_like(Ff))
            Fx = (1. - fac) * F0 + fac * Ff
            dltyp[iel].append(SLOAD)
            dload[iel].append([iedge] + [x for x in Fx])

        # INTERPOLATE DISTRIBUTED LOADS
        for (key, Ff) in self.dloadx.items():
            iel = key
            F0 = self.previous.dloadx.get(key, zeros_like(Ff))
            Fx = (1. - fac) * F0 + fac * Ff
            dltyp[iel].append(DLOAD)
            dload[iel].append(Fx)

        # INTERPOLATE SURFACE FLUXES
        for (key, qf) in self.sfluxx.items():
            iel, iedge = key
            q0 = self.previous.sfluxx.get(key, 0.)
            qn = (1. - fac) * q0 + fac * qf
            dltyp[iel].append(SFLUX)
            dload[iel].append([iedge, qn])

        # INTERPOLATE SURFACE FILMS
        for (key, (Tf, hf)) in self.sfilmx.items():
            iel, iedge = key
            T0, h0 = self.previous.sfilmx.get(key, [0., 0.])
            Too = (1. - fac) * T0 + fac * Tf
            h = (1. - fac) * h0 + fac * hf
            dltyp[iel].append(SFILM)
            dload[iel].append([iedge, Too, h])

        # INTERPOLATE HEAT SOURCES
        for (key, sf) in self.hsrcx.items():
            iel = key
            s0 = self.previous.hsrcx.get(key, zeros_like(sf))
            sx = (1. - fac) * s0 + fac * sf
            dltyp[iel].append(HSRC)
            dload[iel].append(sx)

        return dltyp, dload

    def assign_sload(self, iel, iedge, a):
        self.sloadx[(iel, iedge)] = asarray(a)

    def assign_dload(self, iel, a):
        self.dloadx[iel] = asarray(a)

    def assign_sflux(self, iel, iedge, a):
        self.sfluxx[(iel, iedge)] = asarray(a)

    def assign_sfilm(self, iel, iedge, Too, h):
        self.sfilmx[(iel, iedge)] = [Too, h]

    def assign_hsrc(self, iel, s):
        self.hsrcx[iel] = asarray(s)

    def Frame(self, copy=1):
        if self.frames:
            frame_n = self.frames[-1]
            if not frame_n.converged:
                raise RuntimeError('UNCONVERGED FRAME')
        frame = Frame(self.value)
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = deepcopy(frame_n.field_outputs)
        frame.number = len(self.frames)
        self.frames.append(frame)
        return frame

    @property
    def value(self):
        return self.start + sum([frame.increment for frame in self.frames])

    def copy_from(self, step):
        self.frames[0].field_outputs = deepcopy(step.frames[-1].field_outputs)
        self.dofs[:] = step.dofs
        self.dofx = deepcopy(step.dofx)
        self.cloadx = deepcopy(step.cloadx)
        self.dloadx = deepcopy(step.dloadx)
        self.sloadx = deepcopy(step.sloadx)
        self.sfluxx = deepcopy(step.sfluxx)
        self.sfilmx = deepcopy(step.sfilmx)
        self.hsrcx = deepcopy(step.hsrcx)
        self.predef[:] = step.predef

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

    def RemoveBC(self, nodes, dof):
        self.assign_dof(DIRICHLET, nodes, dof, None)

    def RemoveConcentratedLoad(self, nodes, dof):
        self.assign_dof(NEUMANN, nodes, dof, None)

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

        if amplitude is None:
            # REMOVE THIS BC
            for (i,inode) in enumerate(inodes):
                for j in dofs:
                    I = self.model.dofmap(inode, j)
                    if I is None:
                        logging.warn('INVALID DOF FOR NODE '
                                     '{0}'.format(inode))
                        continue
                    if doftype == DIRICHLET and I in self.dofx:
                        self.dofx.pop(I)
                    elif I in self.cloadx:
                        self.cloadx.pop(I)
            return

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

class Frame(object):
    def __init__(self, start):
        self.start = start
        self._inc = 0.
        self.field_outputs = FieldOutputs()
        self.converged = False

    @property
    def increment(self):
        return self._inc

    @increment.setter
    def increment(self, arg):
        self._inc = arg

    @property
    def value(self):
        return self.start + self.increment

    def FieldOutput(self, ftype, name, position, labels, eleblk=None,
                    ndir=None, nshr=None, ngauss=None,
                    ncomp=None, elements=None, data=None):

        if ftype == SYMTENSOR:
            field = SymmetricTensorField(name, position, labels, ndir, nshr,
                                         eleblk=eleblk, ngauss=ngauss,
                                         elements=elements, data=data)

        elif ftype == VECTOR:
            field = VectorField(name, position, labels, ncomp, eleblk=eleblk,
                                ngauss=ngauss, elements=elements, data=data)

        elif ftype == SCALAR:
            field = ScalarField(name, position, labels, eleblk=eleblk,
                                ngauss=ngauss, elements=elements, data=data)

        elif ftype == TENSOR:
            field = TensorField(name, position, labels, ndir, nshr,
                                eleblk=eleblk, ngauss=ngauss,
                                elements=elements, data=data)

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

    def snapshot(self, dtime, **kwds):

        self.increment = dtime

        # STORE KEYWORDS
        for (kwd, val) in kwds.items():
            if val is None or not val.shape[1]:
                continue
            self.field_outputs[kwd].add_data(val)

        # PASS VALUES FROM SVARS TO THE FRAME OUTPUT
        self.field_outputs.advance()

        self.converged = True

    def foo_bar(self):
        # STORE DEGREES OF FREEDOM
        u, R, temp = self.model.format_dof(dofs)
        if react is not None:
            RF, M, Q = self.model.format_dof(react)

        if temp is not None:
            frame.field_outputs['T'].add_data(temp)
            if react is not None:
                frame.field_outputs['Q'].add_data(Q)
        if u.shape[1]:
            frame.field_outputs['U'].add_data(u)
            if react is not None:
                frame.field_outputs['RF'].add_data(RF)
        if R.shape[1]:
            frame.field_outputs['THETA'].add_data(R)
            if react is not None:
                frame.field_outputs['MOM'].add_data(M)


class SDStep(Step):
    """Base class for stress/displacement steps"""

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
            self.assign_dload(iel, rho*a)

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
        a = asarray(components)
        for iel in ielems:
            self.assign_dload(iel, a)

    def SurfaceLoad(self, surface, components):
        if not is_listlike(components):
            raise UserInputError('EXPECTED SURFACE LOAD VECTOR')
        if len(components) != self.model.dimensions:
            raise UserInputError('EXPECTED {0} SURFACE LOAD '
                                 'COMPONENTS'.format(len(self.model.active_dof)))
        surface = self.model.mesh.find_surface(surface)
        components = [x for x in components]
        for (iel, iedge) in surface:
            self.assign_sload(iel, iedge, components)

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
            components = [x for x in amplitude*n]
            self.assign_sload(iel, iedge, components)

    def Pressure(self, surface, amplitude):
        self.SurfaceLoadN(surface, -amplitude)
