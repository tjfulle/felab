import os
import logging
from numpy import *
import numpy.linalg as la

from .utilities import *
from .constants import *
from .data import StepRepository
from .mesh import Mesh
from .mat import Material
from .exodusii import File

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
    def __init__(self, jobid=None):
        self.jobid = jobid or 'Job-1'
        self.mesh = None
        self.materials = {}
        self.initial_temp = None
        self.fh = None

    @property
    def exofile(self):
        if self.fh is not None:
            return self.fh
        self.fh = File(self.jobid+'.exo', mode='w')
        self.fh.genesis(self.mesh.nodmap, self.mesh.elemap, self.mesh.coord,
                        self.mesh.eleblx, nodesets=self.mesh.nodesets,
                        elemsets=self.mesh.elemsets, sidesets=self.mesh.surfaces)
        return self.fh

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
        if not os.path.isfile(filename):
            raise UserInputError('No such file {0!r}'.format(filename))
        self.mesh = Mesh(filename=filename)
        self._initialize()

    def AbaqusMesh(self, filename):
        """
        Generates a finite element mesh from a Abaqus input file.

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
        if not os.path.isfile(filename):
            raise UserInputError('No such file {0!r}'.format(filename))
        self.mesh = Mesh(filename=filename)
        self._initialize()

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
        if not os.path.isfile(filename):
            raise UserInputError('No such file {0!r}'.format(filename))
        self.mesh = Mesh(filename=filename)
        self._initialize()

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
        self._initialize()

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
        self._initialize()

    def _initialize(self):

        if self.mesh is None:
            raise UserInputError('Mesh must first be created')

        if self.numdim is not None and self.mesh.numdim != self.numdim:
            raise UserInputError('Incorrect mesh dimension')

        self.numele = self.mesh.numele
        self.elements = empty(self.numele, dtype=object)
        self.numnod = self.mesh.numnod
        self.numdim = self.mesh.numdim

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

    @property
    def orphaned_elements(self):
        return [iel for (iel, el) in enumerate(self.elements) if el is None]

    def validate(self, eletyp=None, one=False, fields=None):

        # Validate user input

        if self.orphaned_elements:
            raise UserInputError('All elements must be assigned to an element block')

        if eletyp is not None:
            if not all([isinstance(el, eletyp) for el in self.elements]):
                raise UserInputError('Incorrect element type')

        if one:
            if len(set([type(el) for el in self.elements])) != 1:
                raise UserInputError('Expected only 1 element type')

        # Check consistent loading
        nod_with_bc = {}
        for (i,tag) in enumerate(self.doftags):
            if [t for t in tag[:3] if t == DIRICHLET]:
                nod_with_bc[i] = tag[:3]

        # Now that we know which nodes have bc's, check if any elements have
        # loads applied to that node
        for (iel, el) in enumerate(self.elements):
            if not any(in1d(el.nodes, list(nod_with_bc.keys()))):
                continue
            d = self.sload[iel]
            for (iedge, edge) in enumerate(el.edges):
                if not any(d[iedge]):
                    continue
                # there is a load on this edge
                for nod in el.nodes[edge]:
                    if nod not in nod_with_bc:
                        continue
                    # node is on edge with force, check DOF
                    for (j, v) in enumerate(d[iedge]):
                        if v and nod_with_bc[nod][j] == DIRICHLET:
                            msg = 'Attempting to apply load and displacement '
                            msg += 'on same DOF'
                            logging.warn(msg)
                            #raise UserInputError(msg)

        # Do elements have variable degrees of freedom?
        self.vdof = len(set([el.signature for el in self.elements])) > 1

        if self.vdof:
            # Node freedom association table
            self.nodfat = zeros((self.numnod, MDOF), dtype=int)
            for el in self.elements:
                nfs = el.signature
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
                nfs = eb.eletyp.signature
                for (i, j) in enumerate(nfs):
                    if j:
                        self.active_dof[i] = i
            active_dof = [x for x in self.active_dof if x is not None]

        else:

            # Node freedom association table
            nfs = self.elements[0].signature
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
            self.eftab = array(eftab)

            active_dof = [i for i in range(MDOF) if nfs[i]]

        self.active_dof = array(active_dof)
        self.dofs = zeros((self.numnod, self.active_dof.shape[0]))

        self.validated = True

        node_labels = sorted(self.mesh.nodmap, key=lambda k: self.mesh.nodmap[k])

        self.steps = StepRepository()
        step = self.steps.Step()
        frame = step.frames[0]

        # Node data
        nd = self.numdim
        frame.ScalarField('T', NODE, node_labels)
        frame.VectorField('U', NODE, node_labels, self.numdim)
        if self.initial_temp is not None:
            frame.field_outputs['T'].add_data(self.initial_temp)

        sd = False
        for eb in self.mesh.eleblx:
            if not any(eb.eletyp.signature[:3]):
                continue
            sd = True
            if eb.eletyp.npts is not None:
                args = (INTEGRATION_POINT, eb.labels,
                        eb.eletyp.ndir, eb.eletyp.nshr, eb.eletyp.npts, eb.name)
            else:
                args = (ELEMENT_CENTROID, eb.labels,
                        eb.eletyp.ndir, eb.eletyp.nshr, eb.name)
            frame.SymmetricTensorField('S', *args)
            frame.SymmetricTensorField('E', *args)
            frame.SymmetricTensorField('DE', *args)

        if sd:
            frame.VectorField('R', NODE, node_labels, self.numdim)

    def _element_freedom_table(self):
        eftab = []
        for el in self.elements:
            eft = zeros(el.ndof * el.numnod, dtype=int)
            k, count = 0, 0
            for (i, inode) in enumerate(el.nodes):
                ix, sx = 0, zeros(MDOF, dtype=int)
                nfs1 = self.nodfat[inode]
                for j in range(MDOF):
                    if nfs1[j] > 0:
                        sx[j] = ix
                        ix += 1
                nfs2 = el.signature
                for j in range(MDOF):
                    if nfs2[j] > 0:
                        if nfs1[j] > 0:
                            eft[k] = self.nodfmt[inode] + sx[j]
                            count += 1
                        k += 1
            if all(eft==0.):
                raise UserInputError('Zero entry in eftab for '
                                     'element {0}'.format(el.label))
            eftab.append(eft)
        return eftab

    def request_output_variable(self, name, type, position):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        self.steps[0].FieldOutput(name, type, position)

    def snapshot(self, step=None):
        if step is None:
            step = self.steps.values()[-1]
        if step.written:
            return
        self.exofile.snapshot(step)
        step.written = 1

    def assemble_global_stiffness(self, *args):
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
        ``assemble_global_stiffness`` implements a simplified assembler
        that adopts the following assumptions:

        - nodes are ordered continuously from 0 to :math:`n-1`;
        - there are no multifreedom constraints; and
        - the global stiffness matrix is stored as a full symmetric matrix.

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
        # compute the element stiffness and scatter to global array
        K = zeros((self.numdof, self.numdof))
        for (iel, el) in enumerate(self.elements):
            # Element stiffness
            elearg = [arg[iel] for arg in args]
            Ke = el.stiffness(*elearg)
            eft = self.eftab[iel]
            K[IX(eft, eft)] += Ke
        return K

    def assemble_global_force(self, *args):
        """Assembles the global force matrix

        Parameters
        ----------
        args : tuple
            Arguments to be sent to the individual element force methods.
            Must have length ``numele``

        Returns
        -------
        F : ndarray
            The global force array, with lentgh N, where N is the total number
            of degrees of freedom in the system.

        Notes
        -----
        ``assemble_global_force`` implements a simplified assembler
        that adopts the following assumptions:

        - nodes are ordered continuously from 0 to :math:`n-1`;
        - there are no multifreedom constraints.


        The arguments sent to this function:

        - are sent in by specific application codes;
        - each item in ``args`` must be a list-like variable with length equal to
          the number of elements in the mesh;
        - are passed to individual element force routines as args[i][iel].

        For example, for a plane stress-displacement problem the arguments
        sent are the arrays of distributed and surface loads, ``dloads`` and
        ``sloads``, respectively:

        .. code:: python

           F = self.assemble_global_force(self.dloads, self.sloads)

        This method then calls each element force as:

        .. code:: python

           for (iel, el) in enumerate(self.elements):
              ...
              Fe = el.force(self.dloads[iel], self.sloads[iel])
              ...

        """
        if not self.validated:
            raise PyFEM2Error('Problem must be validated before assembly')

        Q = zeros(self.numdof)
        F = zeros(self.numdof)

        # compute contribution from Neumann boundary
        for n in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[n,j] > 0:
                    if self.doftags[n,j] == NEUMANN:
                        I = self.nodfmt[n] + ix
                        Q[I] = self.dofvals[n,j]
                    ix += 1

        if not args:
            return F, Q

        # compute contribution from element sources and boundary loads
        for (iel, el) in enumerate(self.elements):
            elearg = [arg[iel] for arg in args]
            Fe = el.force(*elearg)
            eft = self.eftab[iel]
            F[eft] += Fe

        return F, Q

    def assemble_global_residual(self, u):
        """
        Assembles the global residual

        """
        # compute contribution from element sources and boundary loads
        R = zeros(self.numdof)
        for (ieb, eb) in self.mesh.eleblx:
            S = self.steps.last.frames[-1].field_outputs[eb.name, 'S']
            for (e, xel) in enumerate(eb.labels):
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                xe = el.xc + u[el.nodes]
                Re = el.residual(xe, S[e])
                R[self.eftab[iel]] += Re
        return R

    def apply_bc(self, K, F, u=None, du=None, fac=1.):
        """
        .. _apply_bc:

        Apply boundary conditions to the global stiffness ``K`` and global
        force ``F``.

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

        if  u is None:  u = zeros(self.numdof)
        if du is None: du = zeros(self.numdof)

        # copy the global arrays
        Kbc, Fbc = K.copy(), F.copy()
        ubc, mask = [], array([False]*self.numdof)

        # Dirichlet boundary conditions
        for i in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[i,j] <= 0:
                    continue
                if self.doftags[i,j] == DIRICHLET:
                    I = self.nodfmt[i] + ix
                    u_cur = u[I] + du[I]
                    ufac = fac * self.dofvals[i,j] - u_cur
                    Fbc -= [K[k,I] * ufac for k in range(self.numdof)]
                    ubc.append(ufac)
                    mask[I] = True
                    Kbc[I,:] = Kbc[:,I] = 0
                    Kbc[I,I] = 1
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
            raise UserInputError('Duplicate material {0!r}'.format(name))
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
                raise UserInputError('Incorrect amplitude length')
            # amplitude is a list of amplitudes
            a = asarray(amplitude)
        dofs = dof if is_listlike(dof) else [dof]
        for (i,inode) in enumerate(inodes):
            for j in dofs:
                if (doftype == DIRICHLET and
                    (self.doftags[inode,j] == NEUMANN and
                     abs(self.doftags[inode,j])>1e-12)):
                    msg = 'Attempting to apply load and displacement '
                    msg += 'on same DOF'
                    raise UserInputError(msg)
                elif doftype == NEUMANN and self.doftags[inode,j]==DIRICHLET:
                    msg = 'Attempting to apply load and displacement '
                    msg += 'on same DOF'
                    raise UserInputError(msg)
                self.doftags[inode, j] = doftype
                self.dofvals[inode, j] = float(a[i])

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def ConcentratedLoad(self, nodes, dof, amplitude=0.):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        self.assign_dof(NEUMANN, nodes, dof, amplitude)

    def GravityLoad(self, region, components):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        orphans = self.orphaned_elements
        if region == ALL:
            ielems = range(self.numele)
        else:
            ielems = [self.mesh.elemap[el] for el in region]
        if not is_listlike(components):
            raise UserInputError('Expected gravity load vector')
        if len(components) != self.numdim:
            raise UserInputError('Expected {0} gravity load '
                                 'components'.format(len(self.active_dof)))
        a = asarray(components)
        for ielem in ielems:
            if ielem in orphans:
                raise UserInputError('Element blocks must be assigned '
                                     'before gravity loads')
            el = self.elements[ielem]
            rho = el.material.density
            if rho is None:
                raise UserInputError('Element material density must be assigned '
                                     'before gravity loads')
            self.dload[ielem] += rho * a

    def DistributedLoad(self, region, components):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        if not is_listlike(components):
            raise UserInputError('Expected distributed load vector')
        if len(components) != self.numdim:
            raise UserInputError('Expected {0} distributed load '
                                 'components'.format(len(self.active_dof)))
        if region == ALL:
            ielems = range(self.numele)
        elif not is_listlike(region):
            ielems = [self.mesh.elemap[region]]
        else:
            ielems = [self.mesh.elemap[el] for el in region]
        if any(in1d(ielems, self.orphaned_elements)):
            raise UserInputError('Element properties must be assigned '
                                 'before DistributedLoad')
        self.dload[ielems] += asarray(components)

    def SurfaceLoad(self, surface, components):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        if not is_listlike(components):
            raise UserInputError('Expected surface load vector')
        if len(components) != self.numdim:
            raise UserInputError('Expected {0} surface load '
                                 'components'.format(len(self.active_dof)))
        surface = self.mesh.find_surface(surface)
        orphans = self.orphaned_elements
        for (iel, iedge) in surface:
            if iel in orphans:
                raise UserInputError('Element properties must be assigned '
                                     'before SurfaceLoad')
            self.sload[iel, iedge, :] = components

    def SurfaceLoadN(self, surface, amplitude):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        surface = self.mesh.find_surface(surface)
        orphans = self.orphaned_elements
        for (iel, edge) in surface:
            # determine the normal to the edge
            if iel in orphans:
                raise UserInputError('Element properties must be assigned '
                                     'before SurfaceLoadN')
            el = self.elements[iel]
            edgenod = el.edges[edge]
            xb = el.xc[edgenod]
            if self.numdim == 2:
                n = normal2d(xb)
            else:
                raise NotImplementedError('3D surface normal')
            self.sload[iel, edge, :] = amplitude * n

    def Pressure(self, surface, amplitude):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        self.SurfaceLoadN(surface, -amplitude)

    # ----------------------------------------------------------------------- #
    # --- HEAT TRANSFER LOADINGS -------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def SurfaceFlux(self, surface, qn):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        surf = self.mesh.find_surface(surface)
        for (iel, edge) in surf:
            self.sflux[iel, edge] = qn

    def SurfaceConvection(self, surface, Too, h):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        surf = self.mesh.find_surface(surface)
        for (iel, edge) in surf:
            self.sfilm[iel, edge, :] = [Too, h]

    def HeatGeneration(self, region, amplitude):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
                raise UserInputError('Heat generation amplitude must have '
                                     'length {0}'.format(len(inodes)))
            a = asarray(amplitude)
        nodmap = dict(zip(inodes, range(inodes.shape[0])))
        for xelem in xelems:
            ielem = self.mesh.elemap[xelem]
            ix = [nodmap[n] for n in self.elements[ielem].nodes]
            self.src[ielem] += a[ix]

    def InitialTemperature(self, nodes, amplitude):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        inodes = self.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, '__call__'):
            # amplitude is a function
            a = amplitude(self.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # create a single amplitude for each node
            a = ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError('Incorrect amplitude length')
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        if elemat not in self.materials:
            raise UserInputError('No such material {0!r}'.format(elemat))
        elemat = self.materials[elemat]
        if blknam.upper() not in self.mesh.element_blocks:
            raise UserInputError('No such element block {0!r}'.format(blknam))
        blk = self.mesh.element_blocks[blknam.upper()]
        blk.eletyp = eletyp
        if eletyp.numnod != blk.elecon.shape[1]:
            raise UserInputError('Node type not consistent with element block')

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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
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
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        self.mesh.ElementSet(name, region)

    def _get_field(self, key):
        key1 = key.lower()
        if key1 in ('u', 'ux', 'uy', 'uz'):
            key1 = 'displ' + key1[1:]
        for (name, field) in self.steps.last.frames[-1].field_outputs.items():
            if key1 == name.lower() or key.lower() == name.lower():
                if field.type != SCALAR:
                    comps = ','.join(key+comp for comp in field.components)
                    msg = 'Non scalar plotting requires components be specified. '
                    msg += 'Try one of {0}'.format(comps)
                    raise UserInputError(msg)
                return field.data
            if key.lower() == name.lower():
                key1 = key.lower()
            if key1 in field.keys:
                if field.position in (ELEMENT, INTEGRATION_POINT):
                    raise NotImplementedError('Plotting element data not done')
                return field.data[:,field.keys.index(key1)]
        raise UserInputError('No such field {0!r}'.format(key))

    def Plot2D(self, deformed=False, color=None, colorby=None, scale=1., **kwds):
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
        if self.numdim != 2:
            raise UserInputError('Plot2D is only applicable to 2D problems')
        xy = array(self.mesh.coord)
        if deformed:
            xy += scale * self.dofs.reshape(xy.shape)
        elecon = []
        for blk in self.mesh.eleblx:
            if (blk.eletyp.numdim, blk.eletyp.numnod) == (2,8):
                raise NotImplementedError('Plotting valid only for linear element')
            else:
                elecon.extend(blk.elecon)

        if colorby is not None and is_stringlike(colorby):
            colorby = self._get_field(colorby)
        return self.mesh.Plot2D(xy=xy, elecon=array(elecon), color=color,
                                colorby=colorby, **kwds)

    def WriteResults(self):
        """Write the finite element results to a file"""
        for (name, step) in self.steps.items():
            self.snapshot(step)
