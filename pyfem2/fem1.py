import os
from numpy import *
import numpy.linalg as la
from utilities import *
from constants import *
from data import TimeStepRepository
from mesh import Mesh
from data import TimeStepRepository
from mat import Material

__all__ = ['FiniteElementModel']

class FiniteElementModel(object):
    """
    The base finite element class

    The ``FiniteElementModel`` object is not intended to be instantiated
    directly. It is a base class from which other finite element models can be
    derived. Specifically, this base class does not provide the ``Solve`` method
    needed to solve a finite element model. The ``Solve`` method must be
    implemented in class' derived from this one.

    """
    active_dof = range(MDOF)
    def __init__(self):
        self.mesh = None
        self.materials = {}
        self.initial_temp = None

    def GenesisMesh(self, filename):
        """
        Generates a finite element mesh from a Genesis file.

        Parameters
        ----------
        filename : str
            The path to a valid Genesis file

        Notes
        -----
        This method calls ``mesh.Mesh`` with the ``filename`` keyword and
        stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

        See Also
        --------
        pyfem2.mesh.Mesh

        """
        assert os.path.isfile(filename), 'No such file {0!r}'.format(filename)
        self.mesh = Mesh(filename=filename)
        self._initialize_geometry()

    def VTKMesh(self, filename):
        """
        Generates a finite element mesh from a vtk .vtu file.

        Parameters
        ----------
        filename : str
            The path to a valid .vtu file

        Notes
        -----
        This method calls ``mesh.Mesh`` with the ``filename`` keyword and
        stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

        See Also
        --------
        pyfem2.mesh.Mesh

        """
        assert os.path.isfile(filename), 'No such file {0!r}'.format(filename)
        self.mesh = Mesh(filename=filename)
        self._initialize_geometry()

    def RectilinearMesh(self, shape, lengths):
        """
        Generates a rectilinear 2D finite element mesh.

        Parameters
        ----------
        shape : tuple
            (nx, ny) where nx is the number elements in :math:`x` and ny
            number of element in :math:`y`.
        lengths : tuple
            (lx, ly) where lx is the length of the mesh in :math:`x` and ny
            is the length of the mesh in :math:`y`.

        Notes
        -----
        This method calls the ``Mesh.RectilinearMesh2D`` class method and
        stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

        See Also
        --------
        pyfem2.mesh.Mesh.RectilinearMesh2D

        """
        self.mesh = Mesh.RectilinearMesh2D(shape, lengths)
        self._initialize_geometry()

    def Mesh(self, **kwds):
        """
        Generates the finite element mesh.

        Notes
        -----
        This method calls the ``mesh.Mesh`` and
        stores the returned mesh as the ``FiniteElementModel.mesh`` attribute.

        See Also
        --------
        pyfem2.mesh.Mesh

        """
        self.mesh = Mesh(**kwds)
        self._initialize_geometry()

    def _initialize_geometry(self):
        assert self.mesh is not None

        if self.numdim is not None and self.mesh.numdim != self.numdim:
            raise ValueError('Incorrect mesh dimension')

        self.numele = self.mesh.numele
        self.elements = empty(self.numele, dtype=object)
        self.numnod = self.mesh.numnod
        self.numdim = self.mesh.numdim

        self.dofs = None
        self.doftags = zeros((self.numnod, MDOF), dtype=int)
        self.dofvals = zeros((self.numnod, MDOF), dtype=float)

        # Containers to hold loads
        self.dload = zeros((self.numele, 2))
        self.sload = zeros((self.numele, self.mesh.maxedge, 2))

        # Heat transfer containers
        #self.src = zeros((self.numele, 3))
        self.src = zeros(self.numele, dtype=object)
        #self.src = zeros(self.numele, self.mesh.maxnod))
        self.sfilm = zeros((self.numele, self.mesh.maxedge, 2))
        self.sflux = zeros((self.numele, self.mesh.maxedge))

        self.validated = False

    def init(self):
        raise NotImplementedError('Application must define init')

    @property
    def orphaned_elements(self):
        return [iel for (iel, el) in enumerate(self.elements) if el is None]

    def validate(self, eletyp=None, one=False):
        # Validate user input
        if self.orphaned_elements:
            raise RuntimeError('All elements must be assigned to an element block')
        if eletyp is not None:
            if not all([isinstance(el, eletyp) for el in self.elements]):
                raise TypeError('Incorrect element type')
        if one:
            if len(set([type(el) for el in self.elements])) != 1:
                raise ValueError('TrussSolution expected only 1 element type')

        # Do elements have variable degrees of freedom?
        self.vdof = len(set([el.signature for el in self.elements])) > 1

        if self.vdof:
            # Node freedom association table
            self.nodfat = zeros((self.numnod, MDOF), dtype=int)
            for el in self.elements:
                nfs = IntegerDigits(el.signature)
                for node in el.nodes:
                    nf = [max(nfs[j], self.nodfat[node,j]) for j in range(MDOF)]
                    self.nodfat[node] = nf

            # Total number of degrees of freedom
            self.numdof = sum(count_digits(p) for p in self.nodfat)

            # Node freedom map table
            self.nodfmt = zeros(self.numnod, dtype=int)
            for i in range(self.numnod-1):
                self.nodfmt[i+1] = self.nodfmt[i] + count_digits(self.nodfat[i])

            # Element freedom table
            self.eftab = self._element_freedom_table()

            # Active degrees of freedom
            self.active_dof = [None] * MDOF
            for eb in self.mesh.eleblx:
                nfs = IntegerDigits(eb.eletyp.signature)
                for (i, j) in enumerate(nfs):
                    if j:
                        self.active_dof[i] = i
            active_dof = [x for x in self.active_dof if x is not None]

        else:

            # Node freedom association table
            nfs = IntegerDigits(self.elements[0].signature)
            self.nodfat = array([nfs for i in range(self.numnod)])
            # Total number of degrees of freedom
            self.numdof = self.numnod * self.elements[0].ndof
            # Node freedom map table
            self.nodfmt = zeros(self.numnod, dtype=int)
            for i in range(self.numnod-1):
                self.nodfmt[i+1] = self.nodfmt[i] + count_digits(self.nodfat[i])
            eftab = []
            for el in self.elements:
                eftab.append(Flatten([[el.ndof*ni+i for i in range(el.ndof)]
                                      for ni in el.nodes]))
            # Element freedom table
            eftab = array(eftab)

            active_dof = [i for i in range(MDOF) if nfs[i]]

        self.steps = TimeStepRepository(self.mesh)
        self.init()
        self.create_step()
        if self.initial_temp is not None:
            self.steps[0].field_outputs['T'].add_data(self.initial_temp)

        self.active_dof = array(active_dof)
        self.dofs = zeros((self.numnod, self.active_dof.shape[0]))

        self.validated = True

    def _element_freedom_table(self):
        eftab = []
        for el in self.elements:
            eft = zeros(el.ndof * el.numnod)
            k, count = 0, 0
            for (i, inode) in enumerate(el.nodes):
                ix, sx = 0, zeros(MDOF, dtype=int)
                nfs1 = self.nodfat[inode]
                for j in range(MDOF):
                    if nfs1[j] > 0:
                        sx[j] = ix
                        ix += 1
                nfs2 = IntegerDigits(el.signature)
                for j in range(MDOF):
                    if nfs2[j] > 0:
                        if nfs1[j] > 0:
                            eft[k] = self.nodfmt[inode] + sx[j]
                            count += 1
                        k += 1
            if eft is None:
                raise ValueError('Zero entry in eftab for '
                                 'element {0}'.format(el.label))
            eftab.append(eft)
        return eftab

    def request_output_variable(self, name, type, position):
        assert self.mesh is not None, 'Mesh must first be created'
        self.steps[0].FieldOutput(name, type, position)

    def create_step(self, dtime=1.):
        assert self.mesh is not None, 'Mesh must first be created'
        # new time step
        ts = self.steps.TimeStep(dtime, copy=1)

    def snapshot(self, **kwds):
        self.steps[-1].add_data(**kwds)

    def _write_steps(self, filename):
        from exodusii import File
        if not filename.endswith(('.e', '.exo')):
            filename += '.exo'
        f = File(filename, mode='w')
        f.genesis(self.mesh.nodmap, self.mesh.elemap, self.mesh.coord,
                  self.mesh.eleblx, nodesets=self.mesh.nodesets,
                  elemsets=self.mesh.elemsets, sidesets=self.mesh.surfaces)
        for step in self.steps:
            f.snapshot(step)

    def assemble_global_stiffness(self, *args):
        """Assembles the global stiffness matrix

        This method is a layer between the other assemblers. It will call the
        appropriate stiffness assembler, depending on whether the finite
        element problem has elements with uniform DOF, or variable DOF.

        Parameters
        ----------
        args : tuple
            Arguments to be sent to individual element stiffness methods.

        Returns
        -------
        K : ndarray
            The NxN global stiffness array, where N is the total number of
            degrees of freedom in the system.

        Notes
        -----
        The arguments sent to this function:

        - are sent in by specific application codes;
        - each item in ``args`` must be a list-like variable with length equal to
          the number of elements in the mesh;
        - are passed to individual element stiffness routines as args[i][iel].

        For example, for a diffusive heat-transfer
        problem the argument sent is the array of film coefficients
        ``sfilm``:

        .. code:: python

           K = self.assemble_global_stiffness(self.sfilm)

        This method then calls each element stiffness as:

        .. code:: python

           for (iel, el) in enumerate(self.elements):
              ...
              Ke = el.stiffness(self.film[iel])
              ...

        """

        if self.vdof:
            return self.assemble_global_stiffness_vdof(*args)
        else:
            return self.assemble_global_stiffness_udof(*args)

    def assemble_global_force(self, *args):
        """Assembles the global force matrix

        This method is a layer between the other assemblers. It will call the
        appropriate force assembler, depending on whether the finite
        element problem has elements with uniform DOF, or variable DOF.

        Parameters
        ----------
        args : tuple
            Arguments to be sent to individual element force methods.

        Returns
        -------
        K : ndarray
            The global force array, with lentgh N, where N is the total number
            of degrees of freedom in the system.

        Notes
        -----
        The arguments sent to this function:

        - are sent in by specific application codes;
        - each item in ``args`` must be a list-like variable with length equal to
          the number of elements in the mesh;
        - are passed to individual element force routines as args[i][iel].

        For example, for a plane stress-displacement problem the arguments
        sent are the arrays of distributed and surface loads, ``dloads`` and
        ``sloads``, respectively:

        .. code:: python

           K = self.assemble_global_force(self.dloads, self.sloads)

        This method then calls each element force as:

        .. code:: python

           for (iel, el) in enumerate(self.elements):
              ...
              Fe = el.force(self.dloads[iel], self.sloads[iel])
              ...

        """
        if self.vdof:
            return self.assemble_global_force_vdof(*args)
        else:
            return self.assemble_global_force_udof(*args)

    def apply_bc(self, K, F):
        """
        .. _apply_bc:

        Apply boundary conditions to the global stiffness ``K`` and global
        force ``F``.

        This method is a layer between the other boundary condition
        application methods. It will call the appropriate method
        depending on whether the finite element problem has elements with
        uniform DOF, or variable DOF.

        Parameters
        ----------
        K : ndarray
            Global stiffness
        F : ndarray
            Global force

        Returns
        -------
        Kbc, Fbc : ndarray
            Boundary condition modified stiffness and force

        Notes
        -----
        Boundary conditions are applied in such a way that ``K`` remains
        symmetric by transferring columns associated with known degrees of
        freedom to ``F``. This method is intended to be called in the
        application code's ``Solve`` method.

        """
        if self.vdof:
            return self.apply_bc_vdof(K, F)
        else:
            return self.apply_bc_udof(K, F)

    def assemble_global_stiffness_udof(self, *args):
        """
        Assembles the global stiffness for problems having uniform DOF.

        A problem with uniform DOF is one in which every element in the mesh
        has the same DOF signature. This method is intended to be called in
        the application code's `Solve` method.

        Parameters
        ----------
        args : tuple
            Arguments to be sent to the individual element stiffness methods.
            Must have length `numele`

        Returns
        -------
        K : ndarray
            The NxN global stiffness array, where N is the total number of degrees
            of freedom in the probem.

        Notes
        -----
        ``assemble_global_stiffness_udof`` implements a simplified assembler
        that adopts the following assumptions:

        - nodes are ordered continuously from 0 to :math:`n-1`;
        - the number and type of degrees of freedom at each node is the same;
        - there are no multifreedom constraints; and
        - the global stiffness matrix is stored as a full symmetric matrix.

        See Also
        --------
        FiniteElementModel.assemble_global_stiffness

        """
        # compute the element stiffness and scatter to global array
        ndof = self.elements[0].ndof
        K = zeros((ndof*self.numnod, ndof*self.numnod))
        for (iel, el) in enumerate(self.elements):
            # Element stiffness
            elearg = [arg[iel] for arg in args]
            Ke = el.stiffness(*elearg)
            # GLOBAL DOF = NUMBER OF DOF PER NODExNODE NUMBER + LOCAL DOF
            eft = Flatten([[ndof*ni+i for i in range(ndof)] for ni in el.nodes])
            neldof = len(eft)
            for i in range(neldof):
                ii = eft[i]
                for j in range(i, neldof):
                    jj = eft[j]
                    K[ii, jj] += Ke[i,j]
                    K[jj, ii] = K[ii, jj]
        return K

    def assemble_global_force_udof(self, *args):
        """
        Assembles the global force for problems having uniform DOF.

        A problem with uniform DOF is one in which every element in the mesh
        has the same DOF signature. This method is intended to be called in
        the application code's ``Solve`` method.

        Parameters
        ----------
        args : tuple
            Arguments to be sent to the individual element force methods.
            Must have length ``numele``

        Returns
        -------
        K : ndarray
            The NxN global stiffness array, where N is the total number of degrees
            of freedom in the probem.

        Notes
        -----
        ``assemble_global_force_udof`` implements a simplified assembler
        that adopts the following assumptions:

        - nodes are ordered continuously from 0 to :math:`n-1`;
        - the number and type of degrees of freedom at each node is the same; and
        - there are no multifreedom constraints.

        See Also
        --------
        FiniteElementModel.assemble_global_force

        """
        doftags = self.doftags[:,self.active_dof]
        dofvals = self.dofvals[:,self.active_dof]
        if doftags.ndim == 1:
            doftags = doftags.reshape(-1,1)
            dofvals = dofvals.reshape(-1,1)
        numnod, ndof = doftags.shape
        # Force contribution on Neummann boundary
        Q = array([dofvals[i,j] if doftags[i,j] == NEUMANN else 0.
                   for i in range(numnod) for j in range(ndof)])
        # compute contribution from element sources and boundary loads
        F = zeros(ndof*numnod)
        for (iel, el) in enumerate(self.elements):
            elearg = [arg[iel] for arg in args]
            Fe = el.force(*elearg)
            # GLOBAL DOF = NUMBER OF DOF PER NODExNODE NUMBER + LOCAL DOF
            eft = Flatten([[ndof*ni+i for i in range(ndof)] for ni in el.nodes])
            for i in range(Fe.shape[0]):
                F[eft[i]] += Fe[i]

        return F, Q

    def apply_bc_udof(self, K, F):
        """
        Applies boundary conditions for problems having uniform DOF.

        A problem with uniform DOF is one in which every element in the mesh
        has the same DOF signature. This method is intended to be called in
        the application code's ``Solve`` method.

        Parameters
        ----------
        K : ndarray
            Global stiffness
        F : ndarray
            Global force

        Returns
        -------
        Kbc, Fbc : ndarray
            Boundary condition modified stiffness and force

        See Also
        --------
        FiniteElementModel.apply_bc

        """
        doftags = self.doftags[:,self.active_dof]
        dofvals = self.dofvals[:,self.active_dof]
        if doftags.ndim == 1:
            doftags = doftags.reshape(-1,1)
            dofvals = dofvals.reshape(-1,1)
        numnod, ndof = doftags.shape
        Kbc, Fbc = K.copy(), F.copy()

        # Dirichlet boundary conditions
        for i in range(numnod):
            for j in range(ndof):
                if doftags[i,j] == DIRICHLET:
                    I = i * ndof + j
                    Fbc -= [K[k,I] * dofvals[i,j] for k in range(numnod*ndof)]
                    Kbc[I, :] = Kbc[:, I] = 0
                    Kbc[I, I] = 1

        # Further modify RHS for Dirichlet boundary
        # This must be done after the loop above.
        for i in range(numnod):
            for j in range(ndof):
                if doftags[i,j] == DIRICHLET:
                    Fbc[i*ndof+j] = dofvals[i,j]

        return Kbc, Fbc

    def assemble_global_stiffness_vdof(self, *args):
        K = zeros((self.numdof, self.numdof))
        for (iel, el) in enumerate(self.elements):
            # Element stiffness
            elearg = [arg[iel] for arg in args]
            Ke = el.stiffness(*elearg)
            eft = self.eftab[iel]
            neldof = len(eft)
            for i in range(neldof):
                ii = eft[i]
                for j in range(i, neldof):
                    jj = eft[j]
                    K[ii,jj] += Ke[i,j]
                    K[jj,ii] = K[ii,jj]
        return K

    def assemble_global_force_vdof(self, *args):
        Q = zeros(self.numdof)
        F = zeros(self.numdof)
        for n in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[n,j] > 0:
                    if self.doftags[n,j] == NEUMANN:
                        ii = self.nodfmt[n] + ix
                        Q[ii] = self.dofvals[n,j]
                    ix += 1
        return F, Q

    def apply_bc_vdof(self, K, F):
        # Apply boundary conditions on copies of K and F
        ubc, mask = [], array([False]*(self.numdof))
        Kbc, Fbc = K.copy(), F.copy()
        for inode in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[inode,j] <= 0:
                    continue
                if self.doftags[inode,j] == DIRICHLET:
                    ii = self.nodfmt[inode] + ix
                    Fbc -= [K[k,ii] * self.dofvals[inode,j]
                            for k in range(self.numdof)]
                    ubc.append(self.dofvals[inode,j])
                    mask[ii] = True
                    Kbc[ii,:] = Kbc[:,ii] = 0
                    Kbc[ii,ii] = 1
                ix += 1
        Fbc[mask] = ubc
        return Kbc, Fbc

    # ----------------------------------------------------------------------- #
    # --- MATERIAL MODELS --------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def Material(self, name):
        """Create an empty material object.

        Parameters
        ----------
        name : str
            The name of the material

        Notes
        -----
        The empty material is put in the self.materials container and can be
        referenced as self.material[name]

        See Also
        --------
        pyfem.mat.Material

        """
        if name in self.materials:
            raise ValueError('Duplicate material {0!r}'.format(name))
        self.materials[name] = Material(name)

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
        assert self.mesh is not None, 'Mesh must first be created'
        self.assign_dof(DIRICHLET, nodes, self.active_dof, 0.)
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
        assert self.mesh is not None, 'Mesh must first be created'
        dof = [x for x in (X,Y,Z) if x in self.active_dof]
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
        assert self.mesh is not None, 'Mesh must first be created'
        self.assign_dof(DIRICHLET, nodes, dof, amplitude)
    PrescribedDOF = PrescribedBC

    def assign_dof(self, doftype, nodes, dof, amplitude):
        inodes = self.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, '__call__'):
            # amplitude is a function
            a = amplitude(self.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # create a single amplitude for each node
            a = ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise ValueError('Incorrect amplitude length')
            # amplitude is a list of amplitudes
            a = asarray(amplitude)
        dofs = dof if is_listlike(dof) else [dof]
        for (i,inode) in enumerate(inodes):
            for j in dofs:
                self.doftags[inode, j] = doftype
                self.dofvals[inode, j] = float(a[i])

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def ConcentratedLoad(self, nodes, dof, amplitude=0.):
        assert self.mesh is not None, 'Mesh must first be created'
        self.assign_dof(NEUMANN, nodes, dof, amplitude)

    def GravityLoad(self, region, components):
        assert self.mesh is not None, 'Mesh must first be created'
        orphans = self.orphaned_elements
        if region == ALL:
            ielems = range(self.numele)
        else:
            ielems = [self.mesh.elemap[el] for el in region]
        if not is_listlike(components):
            raise ValueError('Expected gravity load vector')
        if len(components) != self.numdim:
            raise ValueError('Expected {0} gravity load '
                             'components'.format(len(self.active_dof)))
        a = asarray(components)
        for ielem in ielems:
            if ielem in orphans:
                raise ValueError('Element blocks must be assigned '
                                 'before gravity loads')
            el = self.elements[ielem]
            rho = el.material.density
            if rho is None:
                raise ValueError('Element material density must be assigned '
                                 'before gravity loads')
            self.dload[ielem] += rho * a

    def DistributedLoad(self, region, components):
        assert self.mesh is not None, 'Mesh must first be created'
        if not is_listlike(components):
            raise ValueError('Expected distributed load vector')
        if len(components) != self.numdim:
            raise ValueError('Expected {0} distributed load '
                             'components'.format(len(self.active_dof)))
        if region == ALL:
            ielems = range(self.numele)
        elif not is_listlike(region):
            ielems = [self.mesh.elemap[region]]
        else:
            ielems = [self.mesh.elemap[el] for el in region]
        if any(in1d(ielems, self.orphaned_elements)):
            raise ValueError('Element properties must be assigned '
                             'before DistributedLoad')
        self.dload[ielems] += asarray(components)

    def SurfaceLoad(self, surface, components):
        assert self.mesh is not None, 'Mesh must first be created'
        if not is_listlike(components):
            raise ValueError('Expected surface load vector')
        if len(components) != self.numdim:
            raise ValueError('Expected {0} surface load '
                             'components'.format(len(self.active_dof)))
        surface = self.mesh.find_surface(surface)
        orphans = self.orphaned_elements
        for (iel, iedge) in surface:
            if iel in orphans:
                raise ValueError('Element properties must be assigned '
                                 'before SurfaceLoad')
            self.sload[iel, iedge, :] = components

    def SurfaceLoadN(self, surface, amplitude):
        assert self.mesh is not None, 'Mesh must first be created'
        surface = self.mesh.find_surface(surface)
        orphans = self.orphaned_elements
        for (iel, edge) in surface:
            # determine the normal to the edge
            if iel in orphans:
                raise ValueError('Element properties must be assigned '
                                 'before SurfaceLoadN')
            el = self.elements[iel]
            edgenod = el.edges[edge]
            if self.numdim == 2 and len(edgenod) > 2:
                edgenod = edgenod[[0,-1]]
            xb = el.xc[edgenod]
            if self.numdim == 2:
                n = normal2d(xb)
            else:
                raise NotImplementedError('3D surface normal')
            self.sload[iel, edge, :] = amplitude * n

    def Pressure(self, surface, amplitude):
        assert self.mesh is not None, 'Mesh must first be created'
        self.SurfaceLoadN(surface, -amplitude)

    # ----------------------------------------------------------------------- #
    # --- HEAT TRANSFER LOADINGS -------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def SurfaceFlux(self, surface, qn):
        assert self.mesh is not None, 'Mesh must first be created'
        surf = self.mesh.find_surface(surface)
        for (iel, edge) in surf:
            self.sflux[iel, edge] = qn

    def SurfaceConvection(self, surface, Too, h):
        assert self.mesh is not None, 'Mesh must first be created'
        surf = self.mesh.find_surface(surface)
        for (iel, edge) in surf:
            self.sfilm[iel, edge, :] = [Too, h]

    def HeatGeneration(self, region, amplitude):
        assert self.mesh is not None, 'Mesh must first be created'
        if region == ALL:
            xelems = sorted(self.mesh.elemap.keys())
        else:
            xelems = region
        inodes = []
        for eb in self.mesh.eleblx:
            for (i, xel) in enumerate(eb.labels):
                if xel in xelems:
                    inodes.extend(eb.elecon[i])
        inodes = unique(inodes)
        if hasattr(amplitude, '__call__'):
            # amplitude is a function
            x = self.mesh.coord[inodes]
            a = amplitude(x)
        elif not is_listlike(amplitude):
            a = amplitude * ones(len(inodes))
        else:
            if len(amplitude) != len(inodes):
                raise ValueError('Heat generation amplitude must have '
                                 'length {0}'.format(len(inodes)))
            a = asarray(amplitude)
        nodmap = dict(zip(inodes, range(inodes.shape[0])))
        for xelem in xelems:
            ielem = self.mesh.elemap[xelem]
            ix = [nodmap[n] for n in self.elements[ielem].nodes]
            self.src[ielem] += a[ix]

    def InitialTemperature(self, nodes, amplitude):
        assert self.mesh is not None, 'Mesh must first be created'
        inodes = self.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, '__call__'):
            # amplitude is a function
            a = amplitude(self.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # create a single amplitude for each node
            a = ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise ValueError('Incorrect amplitude length')
            # amplitude is a list of amplitudes
            a = asarray(amplitude)
        self.initial_temp = a

    # ----------------------------------------------------------------------- #
    # --- ELEMENT BLOCKS AND SETS-------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def ElementBlock(self, name, elements):
        """Create an element block and assign elements to it

        Parameters
        ----------
        name : str
            The name of the element block
        elements : int, list, or symbolic constant
            Element label, list of element labels, or symbolic constant

        See Also
        --------
        pyfem2.mesh.Mesh.ElementBlock

        """
        assert self.mesh is not None, 'Mesh must first be created'
        blk = self.mesh.ElementBlock(name, elements)
        return blk

    def AssignProperties(self, blknam, eletyp, elemat, **elefab):
        """Assign properties to elements in an element block

        Parameters
        ----------
        blknam : str
            The name of the element block
        eletyp : object
            The element type (uninstantiated class)
        elemat : str
            The name of the material model
        elefab : dict
            Element fabrication properties

        Notes
        -----
        Before this method is called:

        - elements must be assigned to element blocks; and
        - the material model must be created with the Material method.

        Not all element types require element fabrication properties. For
        example, plane stress-displacement elements require the thickness
        ``t`` be specified but plane heat transfer elements do not

        """
        assert self.mesh is not None, 'Mesh must first be created'
        if elemat not in self.materials:
            raise ValueError('No such material {0!r}'.format(elemat))
        elemat = self.materials[elemat]
        if blknam not in self.mesh.element_blocks:
            raise ValueError('No such elemenb block {0!r}'.format(blknam))
        blk = self.mesh.element_blocks[blknam]
        blk.eletyp = eletyp
        if eletyp.numnod != blk.elecon.shape[1]:
            raise ValueError('Node type not consistent with element block')

        if elefab:
            # element fabrication properties given, make sure there is one
            # property per element
            for (key, val) in elefab.items():
                if not is_listlike(val) or len(val) != len(blk.labels):
                    elefab[key] = [val] * len(blk.labels)
                else:
                    elefab[key] = val
        for (i, xel) in enumerate(blk.labels):
            iel = self.mesh.elemap[xel]
            elenod = blk.elecon[i]
            elecoord = self.mesh.coord[elenod]
            kwds = {}
            for (key, val) in elefab.items():
                kwds[key] = val[i]
            self.elements[iel] = eletyp(xel, elenod, elecoord, elemat, **kwds)

    def NodeSet(self, name, region):
        """Create a node set

        Parameters
        ----------
        name : str
            Name for this element set
        region : int, list, or symbolic constant
            Node label, list of node labels, or symbolic constant

        See Also
        --------
        pyfem2.mesh.Mesh.NodeSet

        """
        assert self.mesh is not None, 'Mesh must first be created'
        self.mesh.NodeSet(name, region)

    def Surface(self, name, surface):
        """Create a surface

        Parameters
        ----------
        name : str
            Name for this element set
        surface : list, list of list, or symbolic constant
            Element/face, list of element/fac, or symbolic constant

        See Also
        --------
        pyfem2.mesh.Mesh.Surface

        """
        assert self.mesh is not None, 'Mesh must first be created'
        self.mesh.Surface(name, surface)

    def ElementSet(self, name, region):
        """Create an element set

        Parameters
        ----------
        name : str
            Name for this element set
        region : int, list, or symbolic constant
            Element label, list of element labels, or symbolic constant

        See Also
        --------
        pyfem2.mesh.Mesh.ElementSet

        """
        assert self.mesh is not None, 'Mesh must first be created'
        self.mesh.ElementSet(name, region)

    def Plot2D(self, deformed=False, color=None, **kwds):
        """Create a 2D plot

        Parameters
        ----------
        deformed : bool, optional {False,True}
            Plot the deformed mesh if True
        color : matplotlib color
        kwds : dict
            kwds passed to pyfem2.mesh.Plot2D

        Returns
        -------
        ax : axes object
            The plot axes

        See Also
        --------
        pyfem2.mesh.Mesh.Plot2D

        """
        assert self.numdim == 2
        xy = array(self.mesh.coord)
        if deformed:
            xy += self.dofs.reshape(xy.shape)
        elecon = []
        for blk in self.mesh.eleblx:
            elecon.extend(blk.elecon)
        return self.mesh.Plot2D(xy=xy, elecon=array(elecon), color=color, **kwds)

    def WriteResults(self, filename):
        """Write the finite element results to a file

        Parameters
        ----------
        filename : str

        Notes
        -----
        If ``filename`` is sent with no extension, the ``.exo`` extension is
        appended.

        """
        self._write_steps(filename)
