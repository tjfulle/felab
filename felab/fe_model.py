import os
import logging
from numpy import *
import numpy.linalg as la

from .utilities import *
from .constants import *
from .mesh import *
from .stage import stage_repository
from .material import Material
from .assembly import vdof, apply_boundary_conditions

__all__ = ['fe_model']

class fe_model(object):
    """The base finite element class"""
    def __init__(self, mesh=None, jobid=None):
        self.jobid = jobid or 'Job-1'
        self.dimensions = None
        self.materials = {}
        self.initial_temp = []
        self.pr_bc = []
        self.fh = None
        self.stages = None
        self._setup = False

        self._mesh = None
        if mesh is not None:
            if not isinstance(mesh, Mesh):
                raise UserInputError('mesh must be a Mesh object')
            self.mesh = mesh

    @property
    def exofile(self):
        if self.fh is not None:
            return self.fh
        self.fh = ExodusFile(self.jobid+'.exo', mode='w')
        self.fh.genesis(self.mesh.nodmap, self.mesh.elemap, self.mesh.coord,
                        self.mesh.element_blocks, nodesets=self.mesh.nodesets,
                        elemsets=self.mesh.elemsets, sidesets=self.mesh.surfaces)
        return self.fh

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):

        if self._mesh is not None:
            logging.warn('MESH ALREADY ASSIGNED, OVERWRITING')

        if not isinstance(mesh, Mesh):
            raise UserInputError('MESH MUST BE A MESH OBJECT')

        self._mesh = mesh
        self.dimensions = self.mesh.dimensions
        self.numele = self.mesh.numele
        self.elements = empty(self.numele, dtype=object)
        self.numnod = self.mesh.numnod
        self._setup = False

    def genesis_mesh(self, filename):
        """
        Generates a finite element mesh from a Genesis file.

        See Also
        --------
        felab.mesh.genesis_mesh

        """
        self.mesh = genesis_mesh(filename=filename)

    def abaqus_mesh(self, filename):
        """
        Generates a finite element mesh from a Abaqus input file.

        See Also
        --------
        felab.mesh.abaqus_mesh

        """
        self.mesh = abaqus_mesh(filename=filename)

    def vtk_mesh(self, filename):
        """
        Generates a finite element mesh from a vtk .vtu file.

        See Also
        --------
        felab.mesh.vtk_mesh

        """
        self.mesh = vtk_mesh(filename=filename)

    def rectilinear_mesh(self, nx=1, ny=1, lx=1, ly=1, shift=None):
        """
        Generates a rectilinear 2D finite element mesh.

        See Also
        --------
        felab.mesh.rectilinear_mesh_2d

        """
        self.mesh = rectilinear_mesh_2d(nx=nx, ny=ny, lx=lx, ly=ly, shift=shift)

    def unit_square_mesh(self, nx=1, ny=1, shift=None):
        """
        Generates a rectilinear 2D finite element mesh.

        See Also
        --------
        felab.mesh.unit_square_mesh

        """
        self.mesh = unit_square_mesh(nx=nx, ny=ny, shift=shift)

    def create_mesh(self, **kwds):
        """
        Generates the finite element mesh.

        See Also
        --------
        felab.mesh.Mesh

        """
        self.mesh = Mesh(**kwds)

    @property
    def dofs(self):
        return self.stages.last.dofs

    @property
    def orphaned_elements(self):
        return [iel for (iel, el) in enumerate(self.elements) if el is None]

    def setup(self):

        if self._setup:
            raise RuntimeError('SETUP MUST BE PERFORMED ONLY ONCE')

        # VALIDATE USER INPUT
        if self.orphaned_elements:
            raise UserInputError('ALL ELEMENTS MUST BE ASSIGNED '
                                 'TO AN ELEMENT BLOCK')

        # CHECK VALIDITY OF ELEMENTS
        self._check_element_validity()

        self.nodfat, self.active_dof = vdof.node_freedom_association_table(
            self.numnod, self.elements, disp=1)

        # TOTAL NUMBER OF DEGREES OF FREEDOM
        self.numdof = vdof.total_degrees_of_freedom(self.nodfat)

        self.nodfmt, self._dofmap = vdof.node_freedom_map_table(self.nodfat, disp=1)

        # ELEMENT FREEDOM TABLE
        self.eftab = vdof.element_freedom_table(self.nodfat, self.nodfmt, self.elements)

        self._setup = True

    def dofmap(self, inode, dof):
        return self._dofmap.get((inode,dof))

    def _check_element_validity(self):
        pass

    def initialize_stages(self):

        node_labels = sorted(self.mesh.nodmap, key=lambda k: self.mesh.nodmap[k])

        self.stages = stage_repository(self)
        stage = self.stages.InitialStage('Stage-0')
        for (nodes, dof) in self.pr_bc:
            stage.assign_prescribed_bc(nodes, dof, amplitude=0.)

        increment = stage.increments[0]

        # NODE DATA
        nd = self.dimensions
        if T in self.active_dof:
            increment.FieldOutput(SCALAR, 'Q', NODE, node_labels)
        increment.FieldOutput(SCALAR, 'T', NODE, node_labels)
        increment.FieldOutput(VECTOR, 'U', NODE, node_labels, ncomp=self.dimensions)
        increment.FieldOutput(VECTOR, 'RF', NODE, node_labels, ncomp=self.dimensions)

        a = in1d((TX,TY,TZ), self.active_dof)
        if any(a):
            n = len([x for x in a if x])
            increment.FieldOutput(VECTOR, 'R', NODE, node_labels, ncomp=n)
            increment.FieldOutput(VECTOR, 'M', NODE, node_labels, ncomp=n)

        if self.initial_temp:
            itemp = self.get_initial_temperature()
            increment.field_outputs['T'].add_data(itemp)

        # ELEMENT DATA
        for eb in self.mesh.element_blocks:

            if not eb.eletyp.variables():
                continue

            ielems = [self.mesh.elemap[xel] for xel in eb.labels]
            elems = self.elements[ielems]

            if eb.eletyp.num_integration():
                position = INTEGRATION_POINT
            else:
                position = ELEMENT_CENTROID

            for variable in eb.eletyp.variables():
                if len(variable) == 2:
                    name, vtype = variable
                    idata = None
                elif len(variable) == 3:
                    name, vtype, idata = variable

                if idata is not None:
                    if idata == 1 and vtype == SYMTENSOR:
                        # IDENTITY
                        idata = array([1.]*elems[0].ndir+[0.]*elems[0].nshr)
                    elif idata == 1 and vtype == TENSOR:
                        idata = eye(elems[0].ndir)

                increment.FieldOutput(vtype, name, position, eb.labels,
                                  ndir=elems[0].ndir, nshr=elems[0].nshr,
                                  eleblk=eb.name, ngauss=elems[0].num_integration(),
                                  elements=elems, ncomp=self.dimensions, data=idata)

        increment.converged = True
        return

    def format_dof(self, dofs):
        # CONSTRUCT DISPLACEMENT AND ROTATION VECTORS
        d1 = len([x for x in self.active_dof if x in (X,Y,Z)])
        u = zeros((self.numnod, d1))

        d2 = len([x for x in self.active_dof if x in (TX,TY,TZ)])
        r = zeros((self.numnod, d2))

        if T in self.active_dof:
            temp = zeros(self.numnod)
        else:
            temp = None

        for n in range(self.numnod):
            ix, ui, ri = 0, 0, 0
            for j in range(MDOF):
                if self.nodfat[n,j] > 0:
                    ii = self.nodfmt[n] + ix
                    if j in (X,Y,Z):
                        u[n,ui] = dofs[ii]
                        ui += 1
                    elif j in (TX,TY,TZ):
                        r[n,ri] = dofs[ii]
                        ri += 1
                    else:
                        temp[n] = dofs[ii]
                    ix += 1
        return u, r, temp

    def snapshot(self, stage=None):

        if stage is None:
            for stage in self.stages.values():
                self.snapshot(stage)
                if stage.written:
                    break
            return

        if stage.written:
            return

        self.exofile.snapshot(stage)
        stage.written = 1

    # ----------------------------------------------------------------------- #
    # --- MATERIAL MODELS --------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def create_material(self, name, **kwargs):
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
        pyfem.material._material.Material

        """
        if name in self.materials:
            raise UserInputError('DUPLICATE MATERIAL {0!r}'.format(name))
        self.materials[name] = Material(name, **kwargs)
        return self.materials[name]

    def assign_prescribed_bc(self, nodes, dof):
        if self.stages is not None:
            raise UserInputError('Boundary conditions must be assigned to stages '
                                 'after creation of first stage')
        self.pr_bc.append((nodes, dof))

    def fix_nodes(self, nodes):
        if self.stages is not None:
            raise UserInputError('Boundary conditions must be assigned to stages '
                                 'after creation of first stage')
        self.pr_bc.append((nodes, ALL))
    fix_dofs = fix_nodes

    def pin_nodes(self, nodes):
        if self.stages is not None:
            raise UserInputError('Boundary conditions must be assigned to stages '
                                 'after creation of first stage')
        self.pr_bc.append((nodes, PIN))

    def assign_initial_temperature(self, nodes, amplitude):
        if self.stages is not None:
            raise UserInputError('Intial temperatures must be assigned '
                                 'before creating first stage')
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
    # --- STAGES ------------------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def unique_stage_name(self):
        i = len(self.stages)
        while 1:
            name = 'Stage-{0}'.format(i)
            if name not in self.stages:
                break
            i += 1
            continue
        return name

    def _validate_stage1(self, nlgeom=False, density=None):
        # VALIDATE INPUT
        for eb in self.mesh.element_blocks:
            iel = self.mesh.elemap[eb.labels[0]]
            el = self.elements[iel]
            if el.material.model.requires:
                if 'nlgeom' in el.material.model.requires and not nlgeom:
                    name = el.material.model.name
                    raise UserInputError('MATERIAL {0!r} REQUIRES '
                                         'nlgeom=True'.format(name.upper()))
            if density and not el.material.density:
                raise UserInputError('STAGE REQUIRES MATERIAL DENSITY')

            if not any(el.signature[0][:3]):
                raise UserInputError('STAGE REQUIRES ELEMENTS WITH '
                                     'DISPLACEMENT DEGREES OF FREEDOM')

    def _validate_stage2(self):
        # VALIDATE INPUT
        for eb in self.mesh.element_blocks:
            iel = self.mesh.elemap[eb.labels[0]]
            el = self.elements[iel]
            if not el.signature[0][T]:
                raise UserInputError('STAGE REQUIRES ELEMENTS WITH '
                                     'TEMPERATURE DEGREE OF FREEDOM')
            if any(el.signature[0][:3]):
                logging.warn('STAGE WILL IGNORE DISPLACEMENT DEGREES OF FREEDOM')

    def create_static_stage(self, name=None, period=1., **kwds):

        if self.stages is None:
            self.setup()
            self.initialize_stages()

        # VALIDATE INPUT
        self._validate_stage1(nlgeom=kwds.get('nlgeom',False))

        if name is None:
            name = self.unique_stage_name()

        if name in self.stages:
            raise UserInputError('Duplicate stage name {0!r}'.format(name))

        stage = self.stages.create_static_stage(name, period, **kwds)
        return stage

    def create_dynamic_stage(self, name=None, period=1., **kwds):

        if period is None:
            raise UserInputError('DYNAMIC STAGE REQUIRES PERIOD')

        if self.stages is None:
            self.setup()
            self.initialize_stages()

        # VALIDATE INPUT
        self._validate_stage1(nlgeom=kwds.get('nlgeom'), density=True)

        if name is None:
            name = self.unique_stage_name()

        if name in self.stages:
            raise UserInputError('Duplicate stage name {0!r}'.format(name))

        stage = self.stages.create_dynamic_stage(name, period, **kwds)
        return stage

    def create_heat_transfer_stage(self, name=None, period=1.):
        if self.stages is None:
            self.setup()
            self.initialize_stages()

        # VALIDATE INPUT
        self._validate_stage2()

        if name is None:
            name = self.unique_stage_name()

        if name in self.stages:
            raise UserInputError('Duplicate stage name {0!r}'.format(name))

        stage = self.stages.create_heat_transfer_stage(name, period=period)
        return stage

    # ----------------------------------------------------------------------- #
    # --- ELEMENT BLOCKS AND SETS-------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def create_element_block(self, name, elements):
        """Create an element block and assign elements to it

        Parameters
        ----------
        name : str
            The name of the element block
        elements : int, list, or symbolic constant
            Element label, list of element labels, or symbolic constant

        See Also
        --------
        felab.mesh.Mesh.element_block

        """
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        blk = self.mesh.create_element_block(name, elements)
        return blk

    def assign_properties(self, blknam, eletyp, elemat, **elefab):
        """Assign properties to elements in an element block

        Parameters
        ----------
        blknam : str
            The name of the element block
        eletyp : object
            The element type (uninstantiated class)
        elemat : str or Material
            The name of the material model, or a material model
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
        if elemat in self.materials:
            elemat = self.materials[elemat]
        elif isinstance(elemat, Material):
            if elemat.name not in self.materials:
                self.materials[elemat.name] = elemat
        else:
            raise UserInputError('NO SUCH MATERIAL {0!r}'.format(elemat))
        for blk in self.mesh.element_blocks:
            if blk.name.upper() == blknam.upper():
                break
        else:
            raise UserInputError('NO SUCH ELEMENT BLOCK {0!r}'.format(blknam))
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

    def create_node_set(self, name, region):
        """Create a node set

        Parameters
        ----------
        name : str
            Name for this element set
        region : int, list, or symbolic constant
            Node label, list of node labels, or symbolic constant

        See Also
        --------
        felab.mesh.Mesh.node_set

        """
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        self.mesh.create_node_set(name, region)

    def create_side_set(self, name, surface):
        """Create a surface

        Parameters
        ----------
        name : str
            Name for this element set
        surface : list, list of list, or symbolic constant
            Element/face, list of element/face, or symbolic constant

        See Also
        --------
        felab.mesh.Mesh.side_set

        """
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        self.mesh.create_side_set(name, surface)

    def create_element_set(self, name, region):
        """Create an element set

        Parameters
        ----------
        name : str
            Name for this element set
        region : int, list, or symbolic constant
            Element label, list of element labels, or symbolic constant

        See Also
        --------
        felab.mesh.Mesh.element_set

        """
        if self.mesh is None:
            raise UserInputError('MESH MUST FIRST BE CREATED')
        self.mesh.create_element_set(name, region)

    def _get_field(self, key):
        key1 = key.lower()
        if key1 in ('u', 'ux', 'uy', 'uz'):
            key1 = 'displ' + key1[1:]
        for (name, field) in self.stages.last.increments[-1].field_outputs.items():
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
            kwds passed to felab.mesh.Plot2D

        Returns
        -------
        ax : axes object
            The plot axes

        See Also
        --------
        felab.mesh.Mesh.Plot2D

        """
        assert self.dimensions == 2
        if self.dimensions != 2:
            raise UserInputError('Plot2D IS ONLY APPLICABLE TO 2D PROBLEMS')
        xy = array(self.mesh.coord)
        if deformed:
            xy += scale * self.stages.last.dofs.reshape(xy.shape)
        elecon = []
        for blk in self.mesh.element_blocks:
            if (blk.eletyp.dimensions, blk.eletyp.nodes) == (2,8):
                raise NotImplementedError('PLOTTING VALID ONLY FOR LINEAR ELEMENT')
            else:
                elecon.extend(blk.elecon)

        if colorby is not None and is_stringlike(colorby):
            colorby = self._get_field(colorby)
        return self.mesh.Plot2D(xy=xy, elecon=array(elecon), color=color,
                                colorby=colorby, **kwds)

    def write_results(self):
        """Write the finite element results to a file"""
        for (name, stage) in self.stages.items():
            self.snapshot(stage)
        self.fh.close()
