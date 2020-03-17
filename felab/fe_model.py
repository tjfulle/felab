import numpy as np

import felab.util.tty as tty
from felab.material import Material
from felab.error import UserInputError
from felab.io import ExodusFile
from felab.mesh import (
    Mesh,
    abaqus_mesh,
    genesis_mesh,
    vtk_mesh,
    unit_square_mesh,
    rectilinear_mesh2d,
)
from felab.step import StepRepository
from felab.util.lang import is_listlike
from felab.assembly import vdof

from felab.constants import (
    X,
    Y,
    Z,
    T,
    TX,
    TY,
    TZ,
    TENSOR,
    INTEGRATION_POINT,
    PIN,
    ELEMENT,
    SCALAR,
    ELEMENT_CENTROID,
    VECTOR,
    NODE,
    SYMTENSOR,
    ALL,
    MDOF,
)

__all__ = ["FEModel"]


class FEModel(object):
    """The base finite element class"""

    def __init__(self, mesh=None, jobid=None):
        self.jobid = jobid or "Job-1"
        self.dimensions = None
        self.materials = {}
        self.initial_temp = []
        self.pr_bc = []
        self.fh = None
        self.steps = None
        self._setup = False

        self._mesh = None
        if mesh is not None:
            if not isinstance(mesh, Mesh):
                raise UserInputError("mesh must be a Mesh object")
            self.mesh = mesh

    @property
    def exofile(self):
        if self.fh is not None:
            return self.fh
        self.fh = ExodusFile(self.jobid + ".exo", mode="w")
        self.fh.genesis(
            self.mesh.nodmap,
            self.mesh.elemap,
            self.mesh.coord,
            self.mesh.element_blocks,
            nodesets=self.mesh.nodesets,
            elemsets=self.mesh.elemsets,
            sidesets=self.mesh.surfaces,
        )
        return self.fh

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):

        if self._mesh is not None:
            tty.warn("MESH ALREADY ASSIGNED, OVERWRITING")

        if not isinstance(mesh, Mesh):
            raise UserInputError("MESH MUST BE A MESH OBJECT")

        self._mesh = mesh
        self.dimensions = self.mesh.dimensions
        self.numele = self.mesh.numele
        self.elements = np.empty(self.numele, dtype=object)
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

    def rectilinear_mesh(
        self, nx=1, ny=1, lx=1, ly=1, shiftx=None, shifty=None, method=None
    ):
        """
        Generates a rectilinear 2D finite element mesh.

        See Also
        --------
        felab.mesh.rectilinear_mesh_2d

        """
        self.mesh = rectilinear_mesh2d(
            nx=nx, ny=ny, lx=lx, ly=ly, shiftx=shiftx, shifty=shifty, method=method
        )

    def unit_square_mesh(self, nx=1, ny=1, shiftx=None, shifty=None, method=None):
        """
        Generates a rectilinear 2D finite element mesh.

        See Also
        --------
        felab.mesh.unit_square_mesh

        """
        self.mesh = unit_square_mesh(
            nx=nx, ny=ny, shiftx=shiftx, shifty=shifty, method=method
        )

    def pt_mesh(self, p, t):
        """
        Generates the finite element mesh.

        See Also
        --------
        felab.mesh.Mesh

        """
        self.mesh = Mesh(p=p, t=t)

    def ne_mesh(self, nodtab, eletab):
        """
        Generates the finite element mesh.

        See Also
        --------
        felab.mesh.Mesh

        """
        self.mesh = Mesh(nodtab=nodtab, eletab=eletab)

    @property
    def dofs(self):
        return self.steps.last.dofs

    @property
    def orphaned_elements(self):
        return [iel for (iel, el) in enumerate(self.elements) if el is None]

    def setup(self):

        if self._setup:
            raise RuntimeError("SETUP MUST BE PERFORMED ONLY ONCE")

        # VALIDATE USER INPUT
        if self.orphaned_elements:
            raise UserInputError("ALL ELEMENTS MUST BE ASSIGNED " "TO AN ELEMENT BLOCK")

        # CHECK VALIDITY OF ELEMENTS
        self._check_element_validity()

        self.nodfat, self.active_dof = vdof.node_freedom_association_table(
            self.numnod, self.elements, disp=1
        )

        # TOTAL NUMBER OF DEGREES OF FREEDOM
        self.numdof = vdof.total_degrees_of_freedom(self.nodfat)

        self.nodfmt, self._dofmap = vdof.node_freedom_map_table(self.nodfat, disp=1)

        # ELEMENT FREEDOM TABLE
        self.eftab = vdof.element_freedom_table(self.nodfat, self.nodfmt, self.elements)

        self._setup = True

    def dofmap(self, inode, dof):
        return self._dofmap.get((inode, dof))

    def _check_element_validity(self):
        pass

    def initialize_steps(self):

        node_labels = sorted(self.mesh.nodmap, key=lambda k: self.mesh.nodmap[k])

        self.steps = StepRepository(self)
        step = self.steps.InitialStep("Step-0")

        # Transfer model boundary conditions to initial step
        for (nodes, dof) in self.pr_bc:
            step.dirichlet_bc(nodes, dof, amplitude=0.0)

        frame = step.frames[0]

        # -- Allocate field outputs for node-based data
        if T in self.active_dof:
            frame.FieldOutput(SCALAR, "Q", NODE, node_labels, description="Heat flux")
        frame.FieldOutput(SCALAR, "T", NODE, node_labels, description="Temperature")
        frame.FieldOutput(
            VECTOR,
            "U",
            NODE,
            node_labels,
            ncomp=self.dimensions,
            description="Displacement",
        )
        frame.FieldOutput(
            VECTOR,
            "RF",
            NODE,
            node_labels,
            ncomp=self.dimensions,
            description="Reaction force",
        )

        # Check for rotation
        a = np.in1d((TX, TY, TZ), self.active_dof)
        if any(a):
            n = len([x for x in a if x])
            frame.FieldOutput(
                VECTOR, "R", NODE, node_labels, ncomp=n, description="Rotation"
            )
            frame.FieldOutput(
                VECTOR, "M", NODE, node_labels, ncomp=n, description="Moment"
            )

        if self.initial_temp:
            itemp = self.get_initial_temperature()
            frame.field_outputs["T"].add_data(itemp)

        # -- Generate the initial element data
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
                        idata = np.array([1.0] * elems[0].ndir + [0.0] * elems[0].nshr)
                    elif idata == 1 and vtype == TENSOR:
                        idata = np.eye(elems[0].ndir)

                frame.FieldOutput(
                    vtype,
                    name,
                    position,
                    eb.labels,
                    ndir=elems[0].ndir,
                    nshr=elems[0].nshr,
                    eleblk=eb.name,
                    ngauss=elems[0].num_integration(),
                    elements=elems,
                    ncomp=self.dimensions,
                    data=idata,
                )

        frame.converged = True
        return

    def format_dof(self, dofs):
        """Construct displacement and rotation vectors from the `dofs`

        `dofs` is a single column vector containing all dofs. This procedure
        separates out components

        """
        d1 = len([x for x in self.active_dof if x in (X, Y, Z)])
        u = np.zeros((self.numnod, d1))

        d2 = len([x for x in self.active_dof if x in (TX, TY, TZ)])
        r = np.zeros((self.numnod, d2))

        if T in self.active_dof:
            temp = np.zeros(self.numnod)
        else:
            temp = None

        for n in range(self.numnod):
            ix, ui, ri = 0, 0, 0
            for j in range(MDOF):
                if self.nodfat[n, j] > 0:
                    ii = self.nodfmt[n] + ix
                    if j in (X, Y, Z):
                        u[n, ui] = dofs[ii]
                        ui += 1
                    elif j in (TX, TY, TZ):
                        r[n, ri] = dofs[ii]
                        ri += 1
                    else:
                        temp[n] = dofs[ii]
                    ix += 1
        return u, r, temp

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

    # ----------------------------------------------------------------------- #
    # --- MATERIAL MODELS --------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def material(self, name, **kwargs):
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
            raise UserInputError("DUPLICATE MATERIAL {0!r}".format(name))
        self.materials[name] = Material(name, **kwargs)
        return self.materials[name]

    def dirichlet_bc(self, nodes, dof):
        if self.steps is not None:
            raise UserInputError(
                "Boundary conditions must be assigned to steps "
                "after creation of first step"
            )
        self.pr_bc.append((nodes, dof))

    def fix_nodes(self, nodes):
        if self.steps is not None:
            raise UserInputError(
                "Boundary conditions must be assigned to steps "
                "after creation of first step"
            )
        self.pr_bc.append((nodes, ALL))

    fix_dofs = fix_nodes

    def pin_nodes(self, nodes):
        if self.steps is not None:
            raise UserInputError(
                "Boundary conditions must be assigned to steps "
                "after creation of first step"
            )
        self.pr_bc.append((nodes, PIN))

    def initial_temperature(self, nodes, amplitude):
        if self.steps is not None:
            raise UserInputError(
                "Intial temperatures must be assigned " "before creating first step"
            )
        self.initial_temp.append((nodes, amplitude))

    def get_initial_temperature(self):
        itemp = np.zeros(self.numnod)
        for (nodes, amplitude) in self.initial_temp:
            inodes = self.mesh.get_internal_node_ids(nodes)
            if hasattr(amplitude, "__call__"):
                # AMPLITUDE IS A FUNCTION
                a = amplitude(self.mesh.coord[inodes])
            elif not is_listlike(amplitude):
                # CREATE A SINGLE AMPLITUDE FOR EACH NODE
                a = np.ones(len(inodes)) * amplitude
            else:
                if len(amplitude) != len(inodes):
                    raise UserInputError("INCORRECT AMPLITUDE LENGTH")
                # AMPLITUDE IS A LIST OF AMPLITUDES
                a = np.asarray(amplitude)
            itemp[inodes] = a
        return itemp

    # ----------------------------------------------------------------------- #
    # --- STEPS ------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def unique_step_name(self):
        i = len(self.steps)
        while 1:
            name = "Step-{0}".format(i)
            if name not in self.steps:
                break
            i += 1
            continue
        return name

    def _validate_step1(self, nlgeom=False, density=None):
        # VALIDATE INPUT
        for eb in self.mesh.element_blocks:
            iel = self.mesh.elemap[eb.labels[0]]
            el = self.elements[iel]
            if el.material.model.requires:
                if "nlgeom" in el.material.model.requires and not nlgeom:
                    name = el.material.model.name
                    raise UserInputError(
                        "MATERIAL {0!r} REQUIRES " "nlgeom=True".format(name.upper())
                    )
            if density and not el.material.density:
                raise UserInputError("STEP REQUIRES MATERIAL DENSITY")

            if not any(el.signature[0][:3]):
                raise UserInputError(
                    "STEP REQUIRES ELEMENTS WITH " "DISPLACEMENT DEGREES OF FREEDOM"
                )

    def _validate_step2(self):
        # VALIDATE INPUT
        for eb in self.mesh.element_blocks:
            iel = self.mesh.elemap[eb.labels[0]]
            el = self.elements[iel]
            if not el.signature[0][T]:
                raise UserInputError(
                    "STEP REQUIRES ELEMENTS WITH " "TEMPERATURE DEGREE OF FREEDOM"
                )
            if any(el.signature[0][:3]):
                tty.warn("STEP WILL IGNORE DISPLACEMENT DEGREES OF FREEDOM")

    def static_step(self, name=None, period=1.0, **kwds):

        if self.steps is None:
            self.setup()
            self.initialize_steps()

        # VALIDATE INPUT
        self._validate_step1(nlgeom=kwds.get("nlgeom", False))

        if name is None:
            name = self.unique_step_name()

        if name in self.steps:
            raise UserInputError("Duplicate step name {0!r}".format(name))

        step = self.steps.static_step(name, period, **kwds)
        return step

    def dynamic_step(self, name=None, period=1.0, **kwds):

        if period is None:
            raise UserInputError("DYNAMIC STEP REQUIRES PERIOD")

        if self.steps is None:
            self.setup()
            self.initialize_steps()

        # VALIDATE INPUT
        self._validate_step1(nlgeom=kwds.get("nlgeom"), density=True)

        if name is None:
            name = self.unique_step_name()

        if name in self.steps:
            raise UserInputError("Duplicate step name {0!r}".format(name))

        step = self.steps.dynamic_step(name, period, **kwds)
        return step

    def heat_transfer_step(self, name=None, period=1.0):
        if self.steps is None:
            self.setup()
            self.initialize_steps()

        # VALIDATE INPUT
        self._validate_step2()

        if name is None:
            name = self.unique_step_name()

        if name in self.steps:
            raise UserInputError("Duplicate step name {0!r}".format(name))

        step = self.steps.heat_transfer_step(name, period=period)
        return step

    # ----------------------------------------------------------------------- #
    # --- ELEMENT BLOCKS AND SETS-------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def element_block(self, *, name, elements):
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
            raise UserInputError("MESH MUST FIRST BE CREATED")
        blk = self.mesh.element_block(name, elements)
        return blk

    def assign_properties(self, *, element_block, element_type, material, **elefab):
        """Assign properties to elements in an element block

        Parameters
        ----------
        element_block : str
            The name of the element block
        element_type : object
            The element type (uninstantiated class)
        material : str or Material
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
            raise UserInputError("MESH MUST FIRST BE CREATED")
        if material in self.materials:
            material = self.materials[material]
        elif isinstance(material, Material):
            if material.name not in self.materials:
                self.materials[material.name] = material
        else:
            raise UserInputError("NO SUCH MATERIAL {0!r}".format(material))
        for blk in self.mesh.element_blocks:
            if blk.name.upper() == element_block.upper():
                break
        else:
            raise UserInputError("NO SUCH ELEMENT BLOCK {0!r}".format(element_block))
        blk.eletyp = element_type
        if element_type.nodes != blk.elecon.shape[1]:
            raise UserInputError("NODE TYPE NOT CONSISTENT WITH ELEMENT BLOCK")

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
            self.elements[iel] = element_type(xel, elenod, elecoord, material, **kwds)

    def node_set(self, name, region):
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
            raise UserInputError("MESH MUST FIRST BE CREATED")
        self.mesh.node_set(name, region)

    def side_set(self, name, surface):
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
            raise UserInputError("MESH MUST FIRST BE CREATED")
        self.side_set(name, surface)

    def element_set(self, name, region):
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
            raise UserInputError("MESH MUST FIRST BE CREATED")
        self.mesh.element_set(name, region)

    def _get_field(self, key):
        key1 = key.lower()
        if key1 in ("u", "ux", "uy", "uz"):
            key1 = "displ" + key1[1:]
        for (name, field) in self.steps.last.frames[-1].field_outputs.items():
            if key1 == name.lower() or key.lower() == name.lower():
                if field.type != SCALAR:
                    comps = ",".join(key + comp for comp in field.components)
                    msg = "NON SCALAR PLOTTING REQUIRES COMPONENTS BE SPECIFIED. "
                    msg += "TRY ONE OF {0}".format(comps)
                    raise UserInputError(msg)
                return field.data
            if key.lower() == name.lower():
                key1 = key.lower()
            if key1 in field.keys:
                if field.position in (ELEMENT, INTEGRATION_POINT):
                    raise NotImplementedError("PLOTTING ELEMENT DATA NOT DONE")
                return field.data[:, field.keys.index(key1)]
        raise UserInputError("NO SUCH FIELD {0!r}".format(key))

    def write_results(self):
        """Write the finite element results to a file"""
        for (name, step) in self.steps.items():
            self.snapshot(step)
        self.fh.close()
