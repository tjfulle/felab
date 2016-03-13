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

def emptywithlists(n):
    a = zeros(n, dtype=object)
    a[:] = [[] for _ in range(n)]
    return a

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
            raise UserInputError('NO SUCH FILE {0!r}'.format(filename))
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
            raise UserInputError('NO SUCH FILE {0!r}'.format(filename))
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
            raise UserInputError('NO SUCH FILE {0!r}'.format(filename))
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
            raise UserInputError('MESH MUST FIRST BE CREATED')

        if self.dimensions is not None and self.mesh.dimensions != self.dimensions:
            raise UserInputError('INCORRECT MESH DIMENSION')

        self.numele = self.mesh.numele
        self.elements = empty(self.numele, dtype=object)
        self.numnod = self.mesh.numnod
        self.dimensions = self.mesh.dimensions

        # CONTAINERS TO HOLD DOF INFORMATION
        self.doftags = zeros((self.numnod, MDOF), dtype=int)
        self.dofvals = zeros((self.numnod, MDOF), dtype=float)

        # CONTAINERS TO HOLD LOADS
        self.dload = emptywithlists(self.numele)
        self.dltyp = emptywithlists(self.numele)

        self._setup = False

    @property
    def orphaned_elements(self):
        return [iel for (iel, el) in enumerate(self.elements) if el is None]

    def setup(self, eletyp=None, one=False, fields=None):

        # VALIDATE USER INPUT

        if self.orphaned_elements:
            raise UserInputError('ALL ELEMENTS MUST BE ASSIGNED '
                                 'TO AN ELEMENT BLOCK')

        if eletyp is not None:
            if not all([isinstance(el, eletyp) for el in self.elements]):
                raise UserInputError('INCORRECT ELEMENT TYPE')

        if one:
            if len(set([type(el) for el in self.elements])) != 1:
                raise UserInputError('EXPECTED ONLY 1 ELEMENT TYPE')

        # CHECK CONSISTENT LOADING
        nod_with_bc = {}
        for (i,tag) in enumerate(self.doftags):
            if [t for t in tag[:3] if t == DIRICHLET]:
                nod_with_bc[i] = tag[:3]

        # NODE FREEDOM ASSOCIATION TABLE
        active_dof = [None] * MDOF
        self.nodfat = zeros((self.numnod, MDOF), dtype=int)
        for el in self.elements:
            for (i, node) in enumerate(el.inodes):
                nfs = el.signature[i]
                nf = [max(nfs[j], self.nodfat[node,j]) for j in range(MDOF)]
                self.nodfat[node] = nf
                for (j, k) in enumerate(nfs):
                    if k:
                        active_dof[j] = j
        active_dof = array([x for x in active_dof if x is not None])

        # TOTAL NUMBER OF DEGREES OF FREEDOM
        self.numdof = sum(count_digits(p) for p in self.nodfat)

        # NODE FREEDOM MAP TABLE
        self.nodfmt = zeros(self.numnod, dtype=int)
        for i in range(self.numnod-1):
            self.nodfmt[i+1] = self.nodfmt[i] + count_digits(self.nodfat[i])

        # ELEMENT FREEDOM TABLE
        self.eftab = self._element_freedom_table()

        self.active_dof = array(active_dof)
        self.dofs = zeros(self.numdof)

        self._setup = True

        node_labels = sorted(self.mesh.nodmap, key=lambda k: self.mesh.nodmap[k])

        # --- ALLOCATE STORAGE FOR SIMULATION DATA
        self.steps = StepRepository()
        step = self.steps.Step()
        frame = step.frames[0]

        # NODE DATA
        nd = self.dimensions
        if 6 in active_dof:
            frame.ScalarField('Q', NODE, node_labels)
        frame.ScalarField('T', NODE, node_labels)
        frame.VectorField('U', NODE, node_labels, self.dimensions)
        frame.VectorField('R', NODE, node_labels, self.dimensions)
        if self.initial_temp is not None:
            frame.field_outputs['T'].add_data(self.initial_temp)

        # ELEMENT DATA
        for eb in self.mesh.eleblx:
            if not eb.eletyp.variables:
                continue
            if eb.eletyp.integration:
                args = (INTEGRATION_POINT, eb.labels,
                        eb.eletyp.ndir, eb.eletyp.nshr,
                        eb.eletyp.integration, eb.name)
            else:
                args = (ELEMENT_CENTROID, eb.labels,
                        eb.eletyp.ndir, eb.eletyp.nshr, eb.name)
            for variable in eb.eletyp.variables:
                frame.SymmetricTensorField(variable, *args)

        frame.converged = True

    def _element_freedom_table(self):
        eftab = []
        for el in self.elements:
            eft = zeros(sum(count_digits(nfs) for nfs in el.signature), dtype=int)
            k, count = 0, 0
            for (i, inode) in enumerate(el.inodes):
                ix, sx = 0, zeros(MDOF, dtype=int)
                nfs1 = self.nodfat[inode]
                for j in range(MDOF):
                    if nfs1[j] > 0:
                        sx[j] = ix
                        ix += 1
                nfs2 = el.signature[i]
                for j in range(MDOF):
                    if nfs2[j] > 0:
                        if nfs1[j] > 0:
                            eft[k] = self.nodfmt[inode] + sx[j]
                            count += 1
                        k += 1
            if all(eft==0.):
                raise UserInputError('ZERO ENTRY IN EFTAB FOR '
                                     'ELEMENT {0}'.format(el.label))
            eftab.append(eft)
        return eftab

    def snapshot(self, step=None):

        if step is None:
            for step in self.steps.values():
                self.snapshot(step)
                if step.written:
                    break
            return

        if step.written:
            return

        self.exofile.snapshot(step)
        step.written = 1

    def assemble(self):
        """
        Assembles the global system of equations

        Returns
        -------
        K : ndarray
            The NxN global stiffness array, where N is the total number of degrees
            of freedom in the probem.

        Notes
        -----
        ``assemble`` implements a simplified assembler that adopts the
        following assumptions:

        - nodes are ordered continuously from 0 to :math:`n-1`;
        - there are no multifreedom constraints; and
        - the global stiffness matrix is stored as a full symmetric matrix.

        flags : ndarray of int
            flags[0] - procedure type
            flags[1] - small or large strain
            flags[2] - general or linear perturbation

        """
        Q = self.external_force_array()

        F = zeros(self.numdof)
        K = zeros((self.numdof, self.numdof))

        # compute the element stiffness and scatter to global array
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            for (e, xel) in enumerate(eb.labels):
                # Element stiffness
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                dload = self.dload[iel]
                dltyp = self.dltyp[iel]
                Ke, Fe = el.response(dltyp, dload)
                eft = self.eftab[iel]
                K[IX(eft, eft)] += Ke
                F[eft] += Fe
                #Ke, Fe = el.response(u[iel], du[iel], time, dtime, fvars[iel],
                #                     istep, iframe, dltyp, dload, ddload, temp,
                #                     flags)

        return K, F, Q

    def external_force_array(self):
        # compute contribution from Neumann boundary
        Q = zeros(self.numdof)
        for n in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[n,j] > 0:
                    if self.doftags[n,j] == NEUMANN:
                        I = self.nodfmt[n] + ix
                        Q[I] = self.dofvals[n,j]
                    ix += 1
        return Q

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
                xe = el.xc + u[el.inodes]
                R[self.eftab[iel]] += el.residual(xe, S[e])
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
            raise UserInputError('DUPLICATE MATERIAL {0!r}'.format(name))
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
                raise UserInputError('INCORRECT AMPLITUDE LENGTH')
            # amplitude is a list of amplitudes
            a = asarray(amplitude)
        dofs = dof if is_listlike(dof) else [dof]
        for (i,inode) in enumerate(inodes):
            for j in dofs:
                if (doftype == DIRICHLET and
                    (self.doftags[inode,j] == NEUMANN and
                     abs(self.doftags[inode,j])>1e-12)):
                    msg = 'ATTEMPTING TO APPLY LOAD AND DISPLACEMENT '
                    msg += 'ON SAME DOF'
                    raise UserInputError(msg)
                elif doftype == NEUMANN and self.doftags[inode,j]==DIRICHLET:
                    msg = 'ATTEMPTING TO APPLY LOAD AND DISPLACEMENT '
                    msg += 'ON SAME DOF'
                    raise UserInputError(msg)
                self.doftags[inode, j] = doftype
                self.dofvals[inode, j] = float(a[i])

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def ConcentratedLoad(self, nodes, dof, amplitude=0.):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        self.assign_dof(NEUMANN, nodes, dof, amplitude)

    def GravityLoad(self, region, components):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        orphans = self.orphaned_elements
        if region == ALL:
            ielems = range(self.numele)
        else:
            ielems = [self.mesh.elemap[el] for el in region]
        if not is_listlike(components):
            raise UserInputError('EXPECTED GRAVITY LOAD VECTOR')
        if len(components) != self.dimensions:
            raise UserInputError('EXPECTED {0} GRAVITY LOAD '
                                 'COMPONENTS'.format(len(self.active_dof)))
        a = asarray(components)
        for iel in ielems:
            if iel in orphans:
                raise UserInputError('ELEMENT BLOCKS MUST BE ASSIGNED '
                                     'BEFORE GRAVITY LOADS')
            el = self.elements[iel]
            rho = el.material.density
            if rho is None:
                raise UserInputError('ELEMENT MATERIAL DENSITY MUST BE ASSIGNED '
                                     'BEFORE GRAVITY LOADS')
            self.dltyp[iel].append(DLOAD)
            self.dload[iel].append(rho*a)

    def DistributedLoad(self, region, components):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        if not is_listlike(components):
            raise UserInputError('EXPECTED DISTRIBUTED LOAD VECTOR')
        if len(components) != self.dimensions:
            raise UserInputError('EXPECTED {0} DISTRIBUTED LOAD '
                                 'COMPONENTS'.format(len(self.active_dof)))
        if region == ALL:
            ielems = range(self.numele)
        elif not is_listlike(region):
            ielems = [self.mesh.elemap[region]]
        else:
            ielems = [self.mesh.elemap[el] for el in region]
        if any(in1d(ielems, self.orphaned_elements)):
            raise UserInputError('ELEMENT PROPERTIES MUST BE ASSIGNED '
                                 'BEFORE DISTRIBUTEDLOAD')
        self.dltyp[ielem].append(DLOAD)
        self.dload[ielem].append(asarray(components))

    def SurfaceLoad(self, surface, components):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        if not is_listlike(components):
            raise UserInputError('EXPECTED SURFACE LOAD VECTOR')
        if len(components) != self.dimensions:
            raise UserInputError('EXPECTED {0} SURFACE LOAD '
                                 'COMPONENTS'.format(len(self.active_dof)))
        surface = self.mesh.find_surface(surface)
        orphans = self.orphaned_elements
        for (iel, iedge) in surface:
            if iel in orphans:
                raise UserInputError('ELEMENT PROPERTIES MUST BE ASSIGNED '
                                     'BEFORE SURFACELOAD')
            self.dltyp[iel].append(SLOAD)
            self.dload[iel].append([iedge]+[x for x in components])

    def SurfaceLoadN(self, surface, amplitude):
        if self.mesh is None:
            raise UserInputError('Mesh must first be created')
        surface = self.mesh.find_surface(surface)
        orphans = self.orphaned_elements
        for (iel, iedge) in surface:
            # determine the normal to the edge
            if iel in orphans:
                raise UserInputError('ELEMENT PROPERTIES MUST BE ASSIGNED '
                                     'BEFORE SURFACELOADN')
            el = self.elements[iel]
            edgenod = el.edges[iedge]
            xb = el.xc[edgenod]
            if self.dimensions == 2:
                n = normal2d(xb)
            else:
                raise NotImplementedError('3D SURFACE NORMAL')
            self.dltyp[iel].append(SLOAD)
            self.dload[iel].append([iedge]+[x for x in amplitude*n])

    def Pressure(self, surface, amplitude):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        self.SurfaceLoadN(surface, -amplitude)

    # ----------------------------------------------------------------------- #
    # --- HEAT TRANSFER LOADINGS -------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def SurfaceFlux(self, surface, qn):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        surf = self.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.dltyp[iel].append(SFLUX)
            self.dload[iel].append([iedge, qn])

    def SurfaceConvection(self, surface, Too, h):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        surf = self.mesh.find_surface(surface)
        for (iel, iedge) in surf:
            self.dltyp[iel].append(SFILM)
            self.dload[iel].append([iedge, Too, h])

    def HeatGeneration(self, region, amplitude):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
                raise UserInputError('HEAT GENERATION AMPLITUDE MUST HAVE '
                                     'LENGTH {0}'.format(len(inodes)))
            a = asarray(amplitude)
        nodmap = dict(zip(inodes, range(inodes.shape[0])))
        for xelem in xelems:
            ielem = self.mesh.elemap[xelem]
            ix = [nodmap[n] for n in self.elements[ielem].inodes]
            self.dltyp[ielem].append(HSRC)
            self.dload[ielem].append(a[ix])

    def InitialTemperature(self, nodes, amplitude):
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        inodes = self.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, '__call__'):
            # amplitude is a function
            a = amplitude(self.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # create a single amplitude for each node
            a = ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError('INCORRECT AMPLITUDE LENGTH')
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
        if elemat not in self.materials:
            raise UserInputError('NO SUCH MATERIAL {0!r}'.format(elemat))
        elemat = self.materials[elemat]
        if blknam.upper() not in self.mesh.element_blocks:
            raise UserInputError('NO SUCH ELEMENT BLOCK {0!r}'.format(blknam))
        blk = self.mesh.element_blocks[blknam.upper()]
        blk.eletyp = eletyp
        if eletyp.nodes != blk.elecon.shape[1]:
            raise UserInputError('NODE TYPE NOT CONSISTENT WITH ELEMENT BLOCK')

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
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
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
            raise UserInputError('MESH MUST FIRST BE CREATED')
        self.mesh.ElementSet(name, region)

    def _get_field(self, key):
        key1 = key.lower()
        if key1 in ('u', 'ux', 'uy', 'uz'):
            key1 = 'displ' + key1[1:]
        for (name, field) in self.steps.last.frames[-1].field_outputs.items():
            if key1 == name.lower() or key.lower() == name.lower():
                if field.type != SCALAR:
                    comps = ','.join(key+comp for comp in field.components)
                    msg = 'NON SCALAR PLOTTING REQUIRES COMPONENTS BE SPECIFIED. '
                    msg += 'TRY ONE OF {0}'.format(comps)
                    raise UserInputError(msg)
                return field.data
            if key.lower() == name.lower():
                key1 = key.lower()
            if key1 in field.keys:
                if field.position in (ELEMENT, INTEGRATION_POINT):
                    raise NotImplementedError('PLOTTING ELEMENT DATA NOT DONE')
                return field.data[:,field.keys.index(key1)]
        raise UserInputError('NO SUCH FIELD {0!r}'.format(key))

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
        assert self.dimensions == 2
        if self.dimensions != 2:
            raise UserInputError('Plot2D IS ONLY APPLICABLE TO 2D PROBLEMS')
        xy = array(self.mesh.coord)
        if deformed:
            xy += scale * self.dofs.reshape(xy.shape)
        elecon = []
        for blk in self.mesh.eleblx:
            if (blk.eletyp.dimensions, blk.eletyp.nodes) == (2,8):
                raise NotImplementedError('PLOTTING VALID ONLY FOR LINEAR ELEMENT')
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
