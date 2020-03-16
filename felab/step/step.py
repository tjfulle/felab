import sys
import numpy as np
import copy as _copy

import felab.util.tty as tty
from felab.error import UserInputError
from felab.util.lang import is_listlike
from felab.util.numeric import emptywithlists, normal2d
from felab.step.data_wharehouse import (
    SymmetricTensorField,
    FieldOutputs,
    ScalarField,
    VectorField,
)
from felab.constants import (
    ALL,
    X,
    Y,
    Z,
    PIN,
    DIRICHLET,
    NEUMANN,
    SCALAR,
    SYMTENSOR,
    VECTOR,
    SFILM,
    SLOAD,
    DLOAD,
    HSRC,
    SFLUX,
)


class load_step(object):
    def __init__(self, model, number, name, previous, period):
        self.model = model
        self.written = 0
        self.ran = False
        self.name = name
        if previous is None:
            self.start = 0.0
            self.value = 0.0
        else:
            self.start = previous.frames[-1].value
            self.value = previous.frames[-1].value
        self.frames = []
        self.Frame(0.0)
        self.period = period
        self.number = number
        self.previous = previous

        self.dofs = np.zeros(self.model.numdof)

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
        self.predef = np.zeros((3, 1, self.model.numnod))

        # --- ALLOCATE STORAGE FOR SIMULATION DATA
        # STATE VARIABLE TABLE
        svtab = []
        nstatev = 0
        for el in self.model.elements:
            if not el.variables():
                svtab.append([])
                continue
            if not el.ndir:
                m = 1
            else:
                m = el.ndir + el.nshr
            m *= len(el.variables())
            if el.num_integration():
                m *= el.num_integration()
            a = [nstatev, nstatev + m]
            svtab.append(slice(*a))
            nstatev += m
        self.svars = np.zeros((2, nstatev))
        self.svtab = svtab
        self._K = None

    def __len__(self):
        return len(self.frames)

    def print_stiffness_structure(
        self, style="ascii", stream=sys.stdout, index_base=0, numeric_fmt="{0: 8f}"
    ):
        if self._K is None:
            return
        m, n = self._K.shape
        assert m == n
        max_int_w = len(str(m))
        if style == "ascii":
            header = " ".join(
                "{0:<{1}d}".format(i + index_base, max_int_w) for i in range(m)
            )
            stream.write(header + "\n")
            for (i, row) in enumerate(self._K):
                for item in row:
                    if abs(item) > 0:
                        stream.write("{0:<{1}s} ".format("X", max_int_w))
                    else:
                        stream.write(" " * (max_int_w + 1))
                stream.write(" {0}\n".format(i + index_base))

        elif style == "latex":
            stream.write("\\setcounter{MaxMatrixCols}{%d}\n" % n)
            stream.write("\\renewcommand{\kbldelim}{[}\n")  # noqa: W605
            stream.write("\\renewcommand{\kbrdelim}{]}\n")  # noqa: W605
            stream.write("\\begin{displaymath}\n  \\kbordermatrix{\n")
            header = "      &" + " & ".join("{0}".format(i) for i in range(m))
            stream.write(header + "\\\\\n")
            for (i, row) in enumerate(self._K):
                s = " & ".join(["\\mathrm{x}" if abs(x) > 0.0 else " " for x in row])
                stream.write("    {0} & {1} \\\\\n".format(i, s))
            stream.write("  }\n\\end{displaymath}\n")

        elif style == "numeric":
            machine_epsilon = np.finfo(float).eps
            fmt = lambda x: numeric_fmt.format(x if abs(x) > machine_epsilon else 0)
            for (i, row) in enumerate(self._K):
                stream.write(" ".join(fmt(x) for x in row))
                stream.write(" {0}\n".format(i + index_base))

        else:
            raise ValueError("Unexpected style {0!r}".format(style))

    @property
    def doftags(self):
        return np.array(sorted(self.dofx), dtype=int)

    def dofvals(self, step_time):

        ix = self.doftags

        # DOFS AT END OF LAST STEP
        X0 = np.array([self.previous.dofx.get(I, 0) for I in ix])

        # DOFS AT END OF THIS STEP
        Xf = np.array([self.dofx[I] for I in ix])

        # INTERPOLATE CONCENTRATED LOAD TO CURRENT TIME
        fac = max(1.0, step_time / self.period)
        return (1.0 - fac) * X0 + fac * Xf

    @property
    def cltags(self):
        return np.array(sorted(self.cloadx), dtype=int)

    def cload(self, step_time):
        # CONCENTRATED LOAD AT END OF LAST STEP
        ix = self.previous.cltags
        Q0 = np.zeros_like(self.dofs)
        Q0[ix] = [self.previous.cloadx[key] for key in ix]

        # CONCENTRATED LOAD AT END OF THIS STEP
        ix = self.cltags
        Qf = np.zeros_like(self.dofs)
        Qf[ix] = [self.cloadx[key] for key in ix]

        # INTERPOLATE CONCENTRATED LOAD TO CURRENT TIME
        fac = max(1.0, step_time / self.period)
        return (1.0 - fac) * Q0 + fac * Qf

    def dload(self, step_time):

        # INTERPOLATES ALL DISTRIBUTED LOADS (BODY AND SURFACE) TO STEP_TIME

        # CONTAINER FOR ALL DLOADS
        dltyp = emptywithlists(self.model.numele)
        dload = emptywithlists(self.model.numele)

        # INTERPOLATION FACTOR
        fac = min(1.0, step_time / self.period)

        # INTERPOLATE SURFACE LOADS
        for (key, Ff) in self.sloadx.items():
            iel, iedge = key
            F0 = self.previous.sloadx.get(key, np.zeros_like(Ff))
            Fx = (1.0 - fac) * F0 + fac * Ff
            dltyp[iel].append(SLOAD)
            dload[iel].append([iedge] + [x for x in Fx])

        # INTERPOLATE DISTRIBUTED LOADS
        for (key, Ff) in self.dloadx.items():
            iel = key
            F0 = self.previous.dloadx.get(key, np.zeros_like(Ff))
            Fx = (1.0 - fac) * F0 + fac * Ff
            dltyp[iel].append(DLOAD)
            dload[iel].append(Fx)

        # INTERPOLATE SURFACE FLUXES
        for (key, qf) in self.sfluxx.items():
            iel, iedge = key
            q0 = self.previous.sfluxx.get(key, 0.0)
            qn = (1.0 - fac) * q0 + fac * qf
            dltyp[iel].append(SFLUX)
            dload[iel].append([iedge, qn])

        # INTERPOLATE SURFACE FILMS
        for (key, (Tf, hf)) in self.sfilmx.items():
            iel, iedge = key
            T0, h0 = self.previous.sfilmx.get(key, [0.0, 0.0])
            Too = (1.0 - fac) * T0 + fac * Tf
            h = (1.0 - fac) * h0 + fac * hf
            dltyp[iel].append(SFILM)
            dload[iel].append([iedge, Too, h])

        # INTERPOLATE HEAT SOURCES
        for (key, sf) in self.hsrcx.items():
            iel = key
            s0 = self.previous.hsrcx.get(key, np.zeros_like(sf))
            sx = (1.0 - fac) * s0 + fac * sf
            dltyp[iel].append(HSRC)
            dload[iel].append(sx)

        return dltyp, dload

    def assign_sload(self, iel, iedge, a):
        self.sloadx[(iel, iedge)] = np.asarray(a)

    def assign_dload(self, iel, a):
        self.dloadx[iel] = np.asarray(a)

    def assign_sflux(self, iel, iedge, a):
        self.sfluxx[(iel, iedge)] = np.asarray(a)

    def assign_sfilm(self, iel, iedge, Too, h):
        self.sfilmx[(iel, iedge)] = [Too, h]

    def assign_hsrc(self, iel, s):
        self.hsrcx[iel] = np.asarray(s)

    def Frame(self, dtime, copy=1):
        frame = Frame(self.value, dtime)
        self.value += dtime
        if self.frames and copy:
            frame_n = self.frames[-1]
            frame.field_outputs = _copy.deepcopy(frame_n.field_outputs)
        frame.number = len(self.frames)
        self.frames.append(frame)
        return frame

    def copy_from(self, step):
        self.frames[0].field_outputs = _copy.deepcopy(
            step.frames[-1].field_outputs
        )
        self.dofs[:] = step.dofs
        self.dofx = _copy.deepcopy(step.dofx)
        self.cloadx = _copy.deepcopy(step.cloadx)
        self.dloadx = _copy.deepcopy(step.dloadx)
        self.sloadx = _copy.deepcopy(step.sloadx)
        self.sfluxx = _copy.deepcopy(step.sfluxx)
        self.sfilmx = _copy.deepcopy(step.sfilmx)
        self.hsrcx = _copy.deepcopy(step.hsrcx)
        self.predef[:] = step.predef
        self.svars[:] = step.svars

    # ----------------------------------------------------------------------- #
    # --- BOUNDARY CONDITIONS ----------------------------------------------- #
    # ----------------------------------------------------------------------- #
    def fix_nodes(self, nodes):
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
        self.assign_dof(DIRICHLET, nodes, ALL, 0.0)

    fix_dofs = fix_nodes

    def remove_bc(self, nodes, dof):
        self.assign_dof(DIRICHLET, nodes, dof, None)

    def remove_concentrated_load(self, nodes, dof):
        self.assign_dof(NEUMANN, nodes, dof, None)

    def assign_prescribed_bc(self, nodes, dof, amplitude=0.0):
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

             self.assign_prescribed_bc(ILO, X, 5.)

        - Assign variable amplitude BC to the :math:`x` displacement of all
          nodes on left side of domain. The variable amplitude function is
          :math:`\Delta_x=y^2`.

          .. code:: python

             fun = lambda x: x[:,1]**2
             self.assign_prescribed_bc(ILO, X, fun)

        """  # noqa: W605
        self.assign_dof(DIRICHLET, nodes, dof, amplitude)

    PrescribedDOF = assign_prescribed_bc

    def assign_dof(self, doftype, nodes, dof, amplitude):

        if dof == ALL:
            dofs = self.model.active_dof
        elif dof == PIN:
            dofs = [x for x in (X, Y, Z) if x in self.model.active_dof]
        elif not is_listlike(dof):
            dofs = [dof]
        else:
            dofs = dof

        inodes = self.model.mesh.get_internal_node_ids(nodes)

        if amplitude is None:
            # REMOVE THIS BC
            for (i, inode) in enumerate(inodes):
                for j in dofs:
                    I = self.model.dofmap(inode, j)  # noqa: E741
                    if I is None:
                        tty.warn("INVALID DOF FOR NODE " "{0}".format(inode))
                        continue
                    if doftype == DIRICHLET and I in self.dofx:
                        self.dofx.pop(I)
                    elif I in self.cloadx:
                        self.cloadx.pop(I)
            return

        if hasattr(amplitude, "__call__"):
            # AMPLITUDE IS A FUNCTION
            a = amplitude(self.model.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # CREATE A SINGLE AMPLITUDE FOR EACH NODE
            a = np.ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError("INCORRECT AMPLITUDE LENGTH")
            # AMPLITUDE IS A LIST OF AMPLITUDES
            a = np.asarray(amplitude)

        for (i, inode) in enumerate(inodes):
            for j in dofs:
                I = self.model.dofmap(inode, j)  # noqa: E741
                if I is None:
                    raise UserInputError("INVALID DOF FOR NODE {0}".format(inode))
                if I in self.cloadx and doftype == DIRICHLET:
                    msg = "ATTEMPTING TO APPLY LOAD AND DISPLACEMENT "
                    msg += "ON SAME DOF"
                    raise UserInputError(msg)
                elif I in self.dofx and doftype == NEUMANN:
                    msg = "ATTEMPTING TO APPLY LOAD AND DISPLACEMENT "
                    msg += "ON SAME DOF"
                    raise UserInputError(msg)
                if doftype == DIRICHLET:
                    self.dofx[I] = float(a[i])
                else:
                    self.cloadx[I] = float(a[i])

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def assign_concentrated_load(self, nodes, dof, amplitude=0.0):
        self.assign_dof(NEUMANN, nodes, dof, amplitude)

    def assign_temperature(self, nodes, amplitude):
        inodes = self.model.mesh.get_internal_node_ids(nodes)
        if hasattr(amplitude, "__call__"):
            # AMPLITUDE IS A FUNCTION
            a = amplitude(self.model.mesh.coord[inodes])
        elif not is_listlike(amplitude):
            # CREATE A SINGLE AMPLITUDE FOR EACH NODE
            a = np.ones(len(inodes)) * amplitude
        else:
            if len(amplitude) != len(inodes):
                raise UserInputError("INCORRECT AMPLITUDE LENGTH")
            # AMPLITUDE IS A LIST OF AMPLITUDES
            a = np.asarray(amplitude)
        self.final_temp = a

    def advance(self, dtime, dofs, react=None, **kwds):
        frame_n = self.frames[-1]
        if not frame_n.converged:
            raise RuntimeError("ATTEMPTING TO UPDATE AN UNCONVERGED FRAME")

        # ADVANCE STATE VARIABLES
        self.svars[0] = self.svars[1]

        # CREATE FRAME TO HOLD RESULTS
        frame = self.Frame(dtime)

        # STORE DEGREES OF FREEDOM
        u, R, temp = self.model.format_dof(dofs)
        if react is not None:
            RF, M, Q = self.model.format_dof(react)

        if temp is not None:
            frame.field_outputs["T"].add_data(temp)
            if react is not None:
                frame.field_outputs["Q"].add_data(Q)
        if u.shape[1]:
            frame.field_outputs["U"].add_data(u)
            if react is not None:
                frame.field_outputs["RF"].add_data(RF)
        if R.shape[1]:
            frame.field_outputs["R"].add_data(R)
            if react is not None:
                frame.field_outputs["M"].add_data(M)

        # STORE KEYWORDS
        for (kwd, val) in kwds.items():
            frame.field_outputs[kwd].add_data(val)

        for (ieb, eb) in enumerate(self.model.mesh.element_blocks):
            if not eb.eletyp.variables():
                continue

            # PASS VALUES FROM SVARS TO THE FRAME OUTPUT
            if eb.eletyp.ndir is not None:
                ntens = eb.eletyp.ndir + eb.eletyp.nshr
            else:
                ntens = None
            m = 1 if not eb.eletyp.num_integration() else eb.eletyp.num_integration()
            n = len(eb.eletyp.variables())
            for (e, xel) in enumerate(eb.labels):
                iel = self.model.mesh.elemap[xel]
                el = self.model.elements[iel]
                if ntens is not None:
                    svars = self.svars[0, self.svtab[iel]].reshape(m, n, ntens)
                else:
                    svars = self.svars[0, self.svtab[iel]].reshape(m, n)
                for (j, variable) in enumerate(el.variables()):
                    name, ftype = variable[:2]
                    frame.field_outputs[eb.name, name].add_data(svars[:, j], e)

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

    def FieldOutput(
        self,
        ftype,
        name,
        position,
        labels,
        eleblk=None,
        ndir=None,
        nshr=None,
        ngauss=None,
        ncomp=None,
        elements=None,
        data=None,
    ):

        if ftype == SYMTENSOR:
            field = SymmetricTensorField(
                name,
                position,
                labels,
                ndir,
                nshr,
                eleblk=eleblk,
                ngauss=ngauss,
                elements=elements,
                data=data,
            )

        elif ftype == VECTOR:
            field = VectorField(
                name,
                position,
                labels,
                ncomp,
                eleblk=eleblk,
                ngauss=ngauss,
                elements=elements,
                data=data,
            )

        elif ftype == SCALAR:
            field = ScalarField(
                name,
                position,
                labels,
                eleblk=eleblk,
                ngauss=ngauss,
                elements=elements,
                data=data,
            )

        if field.eleblk is not None:
            name = (field.eleblk, name)

        self.field_outputs[name] = field

        return field

    def add_data(self, **kwds):
        for (key, value) in kwds.items():
            d = {}
            if isinstance(value, tuple):
                if len(value) == 2:
                    d["ix"] = value[1]
                    value = value[0]
                else:
                    raise ValueError("Unknown add_data option for {0}".format(key))
            self.field_outputs[key].add_data(value, **d)


class sd_step(load_step):
    """Base class for stress/displacement steps"""

    def pin_nodes(self, nodes):
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
        self.assign_dof(DIRICHLET, nodes, PIN, 0.0)

    # ----------------------------------------------------------------------- #
    # --- LOADING CONDITIONS ------------------------------------------------ #
    # ----------------------------------------------------------------------- #
    def assign_gravity_load(self, region, components):
        if region == ALL:
            ielems = range(self.model.numele)
        else:
            ielems = [self.model.mesh.elemap[el] for el in region]
        if not is_listlike(components):
            raise UserInputError("EXPECTED GRAVITY LOAD VECTOR")
        if len(components) != self.model.dimensions:
            raise UserInputError(
                "EXPECTED {0} GRAVITY LOAD "
                "COMPONENTS".format(len(self.model.active_dof))
            )
        a = np.asarray(components)
        for iel in ielems:
            el = self.model.elements[iel]
            rho = el.material.density
            if rho is None:
                raise UserInputError(
                    "ELEMENT MATERIAL DENSITY MUST BE ASSIGNED " "BEFORE GRAVITY LOADS"
                )
            self.assign_dload(iel, rho * a)

    def assign_distributed_load(self, region, components):
        if not is_listlike(components):
            raise UserInputError("EXPECTED DISTRIBUTED LOAD VECTOR")
        if len(components) != self.model.dimensions:
            raise UserInputError(
                "EXPECTED {0} DISTRIBUTED LOAD "
                "COMPONENTS".format(len(self.model.active_dof))
            )
        if region == ALL:
            ielems = range(self.model.numele)
        elif not is_listlike(region):
            ielems = [self.model.mesh.elemap[region]]
        else:
            ielems = [self.model.mesh.elemap[el] for el in region]
        a = np.asarray(components)
        for iel in ielems:
            self.assign_dload(iel, a)

    def assign_surface_load(self, surface, components):
        if not is_listlike(components):
            raise UserInputError("EXPECTED SURFACE LOAD VECTOR")
        if len(components) != self.model.dimensions:
            raise UserInputError(
                "EXPECTED {0} SURFACE LOAD "
                "COMPONENTS".format(len(self.model.active_dof))
            )
        surface = self.model.mesh.find_surface(surface)
        components = [x for x in components]
        for (iel, iedge) in surface:
            self.assign_sload(iel, iedge, components)

    def assign_surface_load_N(self, surface, amplitude):
        surface = self.model.mesh.find_surface(surface)
        for (iel, iedge) in surface:
            # DETERMINE THE NORMAL TO THE EDGE
            el = self.model.elements[iel]
            edgenod = el.edges[iedge]
            xb = el.xc[edgenod]
            if self.model.dimensions == 2:
                n = normal2d(xb)
            else:
                raise NotImplementedError("3D SURFACE NORMAL")
            components = [x for x in amplitude * n]
            self.assign_sload(iel, iedge, components)

    def assign_pressure(self, surface, amplitude):
        self.assign_surface_load_N(surface, -amplitude)
