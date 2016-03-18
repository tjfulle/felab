import os
import logging
from numpy import *
import numpy.linalg as la

from .utilities import *
from .constants import *
from .step_repository import StepRepository
from .mesh.mesh import Mesh
from .mat import Material
from .mesh.exodusii import File

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
    dimensions = None
    def __init__(self, jobid=None):
        self.jobid = jobid or 'Job-1'
        self.mesh = None
        self.materials = {}
        self.initial_temp = []
        self.pr_bc = []
        self.fh = None
        self.steps = None
        self._setup = False

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

    @property
    def dofs(self):
        return self.steps.last.dofs

    def _initialize(self):

        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')

        if self.dimensions is not None and self.mesh.dimensions != self.dimensions:
            raise UserInputError('INCORRECT MESH DIMENSION')

        self.numele = self.mesh.numele
        self.elements = empty(self.numele, dtype=object)
        self.numnod = self.mesh.numnod
        self.dimensions = self.mesh.dimensions

        self._setup = False

    @property
    def orphaned_elements(self):
        return [iel for (iel, el) in enumerate(self.elements) if el is None]

    def setup(self, eletyp=None, one=False):

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
        self.active_dof = array([x for x in active_dof if x is not None])

        # TOTAL NUMBER OF DEGREES OF FREEDOM
        self.numdof = sum(count_digits(p) for p in self.nodfat)

        # NODE FREEDOM MAP TABLE
        self.nodfmt = zeros(self.numnod, dtype=int)
        self._dofmap = {}
        dof = 0
        nodfmt = [0]
        for i in range(self.numnod):
            for (j, k) in enumerate(self.nodfat[i]):
                if not k: continue
                self._dofmap[i,j] = dof
                dof += 1
            nodfmt.append(dof)
        self.nodfmt = array(nodfmt[:-1], dtype=int)

        # ELEMENT FREEDOM TABLE
        self.eftab = self._element_freedom_table()

        self._setup = True

    def dofmap(self, inode, dof):
        return self._dofmap.get((inode,dof))

    def initialize_steps(self):

        node_labels = sorted(self.mesh.nodmap, key=lambda k: self.mesh.nodmap[k])

        self.steps = StepRepository(self)
        step = self.steps.InitialStep('Step-0')
        for (nodes, dof) in self.pr_bc:
            step.PrescribedBC(nodes, dof, amplitude=0.)

        frame = step.frames[0]

        # NODE DATA
        nd = self.dimensions
        if 6 in self.active_dof:
            frame.ScalarField('Q', NODE, node_labels)
        frame.ScalarField('T', NODE, node_labels)
        frame.VectorField('U', NODE, node_labels, self.dimensions)
        frame.VectorField('R', NODE, node_labels, self.dimensions)

        if self.initial_temp:
            itemp = self.get_initial_temperature()
            frame.field_outputs['T'].add_data(itemp)

        # ELEMENT DATA
        for eb in self.mesh.eleblx:
            if not eb.eletyp.variables:
                continue

            ielems = [self.mesh.elemap[xel] for xel in eb.labels]
            elems = self.elements[ielems]

            if eb.eletyp.integration:
                position = INTEGRATION_POINT
            else:
                position = ELEMENT_CENTROID

            for variable in eb.eletyp.variables:
                if elems[0].ndir is None:
                    # SCALAR DATA
                    frame.ScalarField(variable, position, eb.labels,
                                      eleblk=eb.name, ngauss=elems[0].integration,
                                      elements=elems)
                else:
                    # TENSOR DATA
                    frame.SymmetricTensorField(variable, position, eb.labels,
                                               elems[0].ndir, elems[0].nshr,
                                               ngauss=elems[0].integration,
                                               eleblk=eb.name, elements=elems)

        frame.converged = True
        return

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

    def assemble(self, u, du, Q, svtab, svars, dltyp, dload, predef,
                 procedure, step_type, time=array([0.,0.]), dtime=1., period=1.,
                 istep=1, iframe=1, nlgeom=False, ninc=None,
                 cflag=STIFF_AND_FORCE, load_fac=1.):
        """
        Assembles the global system of equations

        Parameters
        ----------
        Returns
        -------
        K : ndarray
            The (N,N) global stiffness array, where N is the total number of degrees
            of freedom in the probem.
        F : ndarray
            The (N,0) global RHS array, where N is the total number of degrees
            of freedom in the probem.

        Notes
        -----
        ``assemble`` implements a simplified assembler that adopts the
        following assumptions:

        - nodes are ordered continuously from 0 to :math:`n-1`;
        - there are no multifreedom constraints; and
        - the global stiffness matrix is stored as a full symmetric matrix.

        """
        procname = get_procname(procedure)
        steptypname = get_steptypname(step_type)
        msg  = 'ASSEMBLING GLOBAL SYSTEM OF EQUATIONS\n      '
        msg += 'PROCEDURE: {0}, STEP TYPE: {1}, NLGEOM: {2}\n      '.format(
            procname, steptypname, nlgeom)
        tf = time[-1] + dtime
        msg += 'STEP: {0}, FRAME: {1}, TIME: {2}'.format(istep, iframe, tf)
        if ninc is not None:
            msg += ', INCREMENT: {0}'.format(ninc)
        logging.debug(msg)

        if cflag not in (STIFF_AND_FORCE, STIFF_ONLY, FORCE_ONLY, LP_OUTPUT):
            raise ValueError('UNKNOWN COMPUTE QUANTITY')

        compute_stiff = cflag in (STIFF_AND_FORCE, STIFF_ONLY)
        compute_force = cflag in (STIFF_AND_FORCE, FORCE_ONLY)
        if compute_stiff:
            K = zeros((self.numdof, self.numdof))

        if compute_force:
            rhs = zeros(self.numdof)

        # INTERPOLATE FIELD VARIABLES
        fac1 = time[1] / (time[0] + period)
        fac2 = (time[1]+dtime) / (time[0] + period)
        x0 = (1. - fac1) * predef[0] + fac1 * predef[1]
        xf = (1. - fac2) * predef[0] + fac2 * predef[1]
        predef_i = array([x0, xf-x0])

        # COMPUTE THE ELEMENT STIFFNESS AND SCATTER TO GLOBAL ARRAY
        for (ieb, eb) in enumerate(self.mesh.eleblx):
            for (e, xel) in enumerate(eb.labels):
                # ELEMENT STIFFNESS
                iel = self.mesh.elemap[xel]
                el = self.elements[iel]
                eft = self.eftab[iel]
                response = el.response(u[eft], du[eft], time, dtime, istep,
                                       iframe, svars[:,svtab[iel]],
                                       dltyp[iel], dload[iel],
                                       predef_i[:,:,el.inodes], procedure, nlgeom,
                                       cflag, step_type, load_fac)

                if cflag == STIFF_AND_FORCE:
                    K[IX(eft, eft)] += response[0]
                    rhs[eft] += response[1]

                elif cflag == STIFF_ONLY:
                    K[IX(eft, eft)] += response

                elif cflag == FORCE_ONLY:
                    rhs[eft] += response

        if cflag == STIFF_AND_FORCE:
            return K, rhs + load_fac * Q

        elif cflag == STIFF_ONLY:
            return K

        elif cflag == FORCE_ONLY:
            return rhs + load_fac * Q

    def apply_bc(self, K, F, doftags, dofvals, u=None, du=None, load_fac=1.):
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

        # COPY THE GLOBAL ARRAYS
        Kbc, Fbc = K.copy(), F.copy()

        # DIRICHLET BOUNDARY CONDITIONS
        for (i, I) in enumerate(doftags):
            u_cur = u[I] + du[I]
            ufac = load_fac * dofvals[i] - u_cur
            Fbc -= [K[k,I] * ufac for k in range(self.numdof)]
            Kbc[I,:] = Kbc[:,I] = 0.
            Kbc[I,I] = 1.
        Fbc[doftags] = dofvals
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

    def PrescribedBC(self, nodes, dof):
        if self.steps is not None:
            raise UserInputError('Boundary conditions must be assigned to steps '
                                 'after creation of first step')
        self.pr_bc.append((nodes, dof))

    def FixNodes(self, nodes):
        if self.steps is not None:
            raise UserInputError('Boundary conditions must be assigned to steps '
                                 'after creation of first step')
        self.pr_bc.append((nodes, ALL))
    FixDOF = FixNodes

    def PinNodes(self, nodes):
        if self.steps is not None:
            raise UserInputError('Boundary conditions must be assigned to steps '
                                 'after creation of first step')
        self.pr_bc.append((nodes, PIN))

    def InitialTemperature(self, nodes, amplitude):
        if self.steps is not None:
            raise UserInputError('Intial temperatures must be assigned '
                                 'before creating first step')
        self.initial_temp.append((nodes, amplitude))

    def get_initial_temperature(self):
        itemp = zeros(self.numnod)
        for (nodes, amplitude) in self.initial_temp:
            inodes = self.mesh.get_internal_node_ids(nodes)
            if hasattr(amplitude, '__call__'):
                # AMPLITUDE IS A FUNCTION
                a = amplitude(self.mesh.coord[inodes])
            elif not is_listlike(amplitude):
                # CREATE A SINGLE AMPLITUDE FOR EACH NODE
                a = ones(len(inodes)) * amplitude
            else:
                if len(amplitude) != len(inodes):
                    raise UserInputError('INCORRECT AMPLITUDE LENGTH')
                # AMPLITUDE IS A LIST OF AMPLITUDES
                a = asarray(amplitude)
            itemp[inodes] = a
        return itemp

    # ----------------------------------------------------------------------- #
    # --- STEPS ------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def unique_step_name(self):
        i = len(self.steps)
        while 1:
            name = 'Step-{0}'.format(i)
            if name not in self.steps:
                break
            i += 1
            continue
        return name

    def StaticStep(self, name=None, period=1., increments=None, maxiters=10,
                   nlgeom=False, solver=None):
        if self.steps is None:
            self.setup()
            self.initialize_steps()
        if name is None:
            name = self.unique_step_name()
        if name in self.steps:
            raise UserInputError('Duplicate step name {0!r}'.format(name))
        step = self.steps.StaticStep(name, period, increments,
                                     maxiters, nlgeom, solver)
        return step

    def HeatTransferStep(self, name=None, period=1.):
        if self.steps is None:
            self.setup()
            self.initialize_steps()
        if name is None:
            name = self.unique_step_name()
        if name in self.steps:
            raise UserInputError('Duplicate step name {0!r}'.format(name))
        step = self.steps.HeatTransferStep(name, period=period)
        return step

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
            # ELEMENT FABRICATION PROPERTIES GIVEN, MAKE SURE THERE IS ONE
            # PROPERTY PER ELEMENT
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
            Element/face, list of element/face, or symbolic constant

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
            xy += scale * self.steps.last.dofs.reshape(xy.shape)
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
        self.fh.close()