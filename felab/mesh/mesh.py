import numpy as np

from felab.constants import ALL, ILO, IHI, JLO, JHI, KLO, KHI, BOUNDARY

from felab.error import UserInputError
from felab.util.lang import is_listlike, is_stringlike
from felab.elemlib import element_family

from felab.mesh import aba
from felab.mesh import vtk
from felab.mesh import exodusii
from felab.mesh.element_block import element_block


__all__ = ["Mesh"]


class Mesh(object):
    """Class for creating a finite element mesh.

    Parameters
    ----------
    nodtab : list of list
        ``nodtab[n]`` is a list describing the nth node where
        ``nodtab[n][0]`` is the node label of the nth node and
        ``nodtab[n][1:]`` are the nodal coordinates of the nth node
    eletab : list of list
        ``eletab[e]`` is a list describing the eth element where
        ``eletab[e][0]`` is the element label of the eth element and
        ``eletab[e][1:]`` are the nodes forming the eth element, the nodal
        labels must correpond to the node labels in ``nodtab``
    filename : str
        Name of a genesis or vtu file containing the mesh
    p : ndarray of float
        Array of nodal coordinates
    t : ndarray of int
        Array of element connectivity

    Notes
    -----
    The nodtab/eletab, filename, and p/t arguments are mutually exclusive.

    If ``p`` and ``t`` are specified, node labels are assigned to nodes in the
    order they appear in the nodal coordinates array ``p`` and take values of
    ``0`` to ``n-1``, where ``n`` is the length of ``p``. The element
    connectivity array ``t`` must use this numbering scheme. This method is
    useful for creating meshes using triangulations created by the
    ``distmesh2d`` module.

    Examples
    --------

    In the following examples, a simple mesh is created of a unit square
    centered at the origin with node and element numbers as follows:

    .. image:: unitsquare1.png
       :scale: 35
       :align: center

    Method 1: specifying ``nodtab`` and ``eletab``:

    >>> nodtab = [[10, -.5, -.5], [20, .5, -.5], [30, .5, .5], [40, -.5, .5]]
    >>> eletab = [[100, 10, 20, 30, 40]]
    >>> mesh = Mesh(nodtab=nodtab, eletab=eletab)

    Method 2: specifying ``filename``

    >>> mesh = Mesh(filename='unitsquare.g')

    Method 3: specifying ``p`` and ``t``:

    >>> p = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
    >>> t = np.array([[0, 1, 2, 3]])
    >>> mesh = Mesh(p=p, t=t)

    """

    def __init__(self, nodtab=None, eletab=None, filename=None, p=None, t=None):
        o1 = nodtab is not None and eletab is not None
        o2 = filename is not None
        o3 = p is not None and t is not None
        if not (o1 or o2 or o3):
            raise UserInputError("nodtab/eletab, filename, or p/t required")
        elif len([x for x in (o1, o2, o3) if x]) > 1:
            raise UserInputError(
                "nodtab/eletab, filename, an p/t are " "mutually exclusive"
            )
        if o1:
            # Format nodes and elements
            data = self.parse_nod_and_elem_tables(nodtab, eletab)
            self.init1(*data)

        elif o2:
            self.init_from_file(filename)

        elif o3:
            p, t = np.asarray(p, dtype=float), np.asarray(t, dtype=int)
            nodmap = dict(zip(range(p.shape[0]), range(p.shape[0])))
            eletab = dict([(i, t[i]) for i in range(t.shape[0])])
            self.init1(nodmap, p, eletab)

        else:
            raise UserInputError("nodtab/eletab or filename required")

        self._bndry_nod2 = None
        self._free_edges = None

    @property
    def unassigned(self):
        return self.num_assigned != self.numele

    def connectivity(self):
        M = len(self.eletab)
        N = max([len(x) for x in self.eletab.values()])
        elecon = np.zeros((M, N), dtype=int)

        a = sorted(self.eletab.keys())
        for (i, e) in enumerate(a):
            ec = self.eletab[e]
            elecon[i, : len(ec)] = ec
        return elecon

    def init1(self, nodmap, coord, eletab, nodesets=None, elemsets=None, surfaces=None):
        self.nodmap = nodmap
        self.coord = np.asarray(coord)
        self.numnod, self.dimensions = self.coord.shape

        self.eletab = eletab
        self.numele = len(self.eletab)
        self.elemap = dict(zip(self.eletab.keys(), [None] * self.numele))
        self.ielemap = {}

        self.nodesets = nodesets or {}
        self.elemsets = elemsets or {}
        self.surfaces = surfaces or {}

        self.element_blocks = []

        # maximum number of edges on any one element
        self.maxedge = max(
            [
                len(element_family(self.dimensions, len(nodes)).edges)
                for nodes in self.eletab.values()
            ]
        )

        # Number of elements assigned to a block
        self.num_assigned = 0

    def init2(
        self, nodmap, coord, elemap, element_blocks, nodesets, elemsets, sidesets
    ):
        self.nodmap = nodmap
        self.coord = np.asarray(coord)
        self.numnod, self.dimensions = self.coord.shape

        self.elemap = elemap
        self.numele = len(self.elemap)
        self.ielemap = np.array(sorted(elemap.keys(), key=lambda k: elemap[k]))

        self.nodesets, self.elemsets, self.surfaces = {}, {}, {}
        for (name, nodes) in nodesets.items():
            self.nodesets[name.upper()] = nodes
        for (name, elems) in elemsets.items():
            self.elemsets[name.upper()] = elems
        for (name, surf) in sidesets.items():
            self.surfaces[name.upper()] = surf

        self.element_blocks = element_blocks
        self.maxedge = max([len(eb.elefam.edges) for eb in element_blocks])
        self.num_assigned = self.numele

    def init_from_file(self, filename):
        if filename.endswith((".vtk", ".vtu")):
            data = self.parse_nod_and_elem_tables(*vtk.read_mesh(filename))
            self.init1(*data)

        elif filename.endswith((".exo", ".g", ".gen")):
            exo = exodusii.File(filename, mode="r")
            self.init2(
                exo.nodmap,
                exo.coord,
                exo.elemap,
                exo.element_blocks,
                exo.nodesets,
                exo.elemsets,
                exo.sidesets,
            )

        elif filename.endswith(".inp"):
            data = aba.ReadInput(filename)
            nodtab, eletab, nodesets, elemsets, surfaces, element_blocks = data
            nodmap, coord, eletab1 = self.parse_nod_and_elem_tables(nodtab, eletab)
            self.init1(nodmap, coord, eletab1)
            for (name, (et, elems)) in element_blocks.items():
                self.element_block(name, elems)
            for (name, nodes) in nodesets.items():
                self.node_set(name, nodes)
            for (name, elems) in elemsets.items():
                self.element_set(name, elems)
            for (name, surf) in surfaces.items():
                self.side_set(name, surf)

        else:
            raise UserInputError("Unknown file type")

    def get_internal_node_ids(self, label):
        if is_stringlike(label) and label.upper() in self.nodesets:
            return self.nodesets[label.upper()]
        elif isinstance(label, int):
            inodes = [self.nodmap[label]]
        elif label == ALL:
            inodes = range(self.numnod)
        elif label == BOUNDARY:
            inodes = self.boundary_nodes()
        elif label in (ILO, IHI, JLO, JHI, KLO, KHI):
            inodes = self.nodes_in_rectilinear_region(label)
        elif is_stringlike(label):
            raise UserInputError("No such node set {0!r}".format(label))
        else:
            inodes = [self.nodmap[xn] for xn in label]
        return np.array(inodes, dtype=int)

    @staticmethod
    def parse_nod_and_elem_tables(nodtab, eletab):
        """

        .. _parse_nodes_and_elements:

        Format the node map, coordinates array, element map, types, and
        connectivity table

        Parameters
        ----------
        nodtab : list of list
            ``nodtab[n]`` is a list describing the nth node where
            ``nodtab[n][0]`` is the node label of the nth node and
            ``nodtab[n][1:]`` are the nodal coordinates of the nth node
        eletab : list of list
            ``eletab[e]`` is a list describing the eth element where
            ``eletab[e][0]`` is the element label of the eth element and
            ``eletab[e][1:]`` are the nodes forming the eth element, the
            nodal labels must correpond to the node labels in ``nodtab``

        Returns
        -------
        nodmap : dict
            nodmap[xn] is the internal node number of node labeled xn
        coord : ndarray
            Nodal coordinates.  coord[:,0] are the x, coord[:,1] the y, etc.
        eletab1 : dict
            eletab1[xel] are the internal node numbers of element labeled xel

        """
        # Basic node info
        nodmap = {}
        numnod, maxdim = len(nodtab), 0
        for (inode, noddef) in enumerate(nodtab):
            if noddef[0] in nodmap:
                raise UserInputError("Duplicate node label: {0}".format(noddef[0]))
            nodmap[noddef[0]] = inode
        maxdim = max(maxdim, len(noddef[1:]))
        coord = np.zeros((numnod, maxdim))
        for (inode, noddef) in enumerate(nodtab):
            coord[inode, : len(noddef[1:])] = noddef[1:]

        # Basic element info
        eletab1 = {}
        for (iel, eledef) in enumerate(eletab):
            if eledef[0] in eletab1:
                raise UserInputError(
                    "Duplicate element label: " "{0}".format(eledef[0])
                )
            eletab1[eledef[0]] = [nodmap[n] for n in eledef[1:]]

        return nodmap, coord, eletab1

    def find_edges(self):
        # Find edges
        if self.unassigned:
            raise UserInputError(
                "Mesh element operations require all elements be "
                "assigned to an element block"
            )
        if self._free_edges is not None:
            return self._free_edges

        edges = {}
        for (ieb, eb) in enumerate(self.element_blocks):
            for (j, elenod) in enumerate(eb.elecon):
                iel = self.elemap[eb.labels[j]]
                for (i, ix) in enumerate(eb.elefam.edges):
                    edge = tuple(elenod[ix])
                    edges.setdefault(tuple(sorted(edge)), []).append((iel, i) + edge)
        self._free_edges = edges.values()
        return self._free_edges

    def boundary_nodes2(self):
        if self._bndry_nod2 is not None:
            return self._bndry_nod2
        edges = self.find_edges()
        nodes = []
        for edge in edges:
            if len(edge) == 1:
                # edge not shared by multiple elements
                nodes.extend(edge[0][2:])
        self._bndry_nod2 = np.unique(sorted(nodes))
        return self._bndry_nod2

    def boundary_nodes(self):
        if self.dimensions == 1:
            return [np.argmin(self.coord), np.argmax(self.coord)]
        elif self.dimensions == 2:
            return self.boundary_nodes2()
        raise UserInputError("3D meshes not supported")

    def nodes_in_rectilinear_region(self, region, tol=1e-6):
        # Find nodes in region
        if region not in (ILO, IHI, JLO, JHI, KLO, KHI):
            raise UserInputError("unknown region {0!r}".format(region))
        axis, fun = {
            ILO: (0, np.amin),
            IHI: (0, np.amax),
            JLO: (1, np.amin),
            JHI: (1, np.amax),
            KLO: (2, np.amin),
            KHI: (2, np.amax),
        }[region]
        xpos = fun(self.coord[:, axis])
        return np.where(abs(self.coord[:, axis] - xpos) < tol)[0]

    def find_surface(self, region):
        if self.unassigned:
            raise UserInputError(
                "Mesh element operations require all elements be "
                "assigned to an element block"
            )
        if is_stringlike(region) and region.upper() in self.surfaces:
            return self.surfaces[region.upper()]
        if region in (ILO, IHI, JLO, JHI, KLO, KHI):
            if self.dimensions == 1:
                return self.find_surface1(region)
            if self.dimensions == 2:
                return self.find_surface2(region)
        elif is_stringlike(region):
            raise UserInputError("No such surface {0!r}".format(region))
        # Check if region is a surface
        if not is_listlike(region):
            raise UserInputError("Unrecognized surface: {0}".format(region))
        surface = []
        region = np.asarray(region)
        if region.ndim == 1:
            region = region[np.newaxis, :]
        for item in region:
            elem, edge = [int(x) for x in item]
            if is_listlike(elem):
                surface.extend([(self.elemap[e], edge) for e in elem])
            else:
                surface.append((self.elemap[elem], edge))
        return surface

    def find_surface1(self, region):
        if region not in (ILO, IHI):
            raise UserInputError("Incorrect 1D region")
        if region == ILO:
            node = np.argmin(self.coord)
        else:
            node = np.argmax(self.coord)
        for (ieb, eb) in enumerate(self.element_blocks):
            for (e, elenod) in enumerate(eb.elecon):
                if node in elenod:
                    i = 0 if node == elenod[0] else 1
                    return [self.elemap[eb.labels[e]], i]

    def find_surface2(self, region):
        nodes = self.nodes_in_rectilinear_region(region)
        surface = []
        for (ieb, eb) in enumerate(self.element_blocks):
            for (e, elenod) in enumerate(eb.elecon):
                w = np.where(np.in1d(elenod, nodes))[0]
                if len(w) < 2:
                    continue
                w = tuple(sorted(w))
                for (ie, edge) in enumerate(eb.elefam.edges):
                    if tuple(sorted(edge)) == w:
                        # the internal node numbers match, this is the edge
                        surface.append((self.elemap[eb.labels[e]], ie))
        return np.array(surface)

        for (e, c) in enumerate(self.elecon):
            c = c[: self.elefam[e].nodes]
            w = np.where(np.in1d(c, nodes))[0]
            if len(w) < 2:
                continue
            w = tuple(sorted(w))
            for (ie, edge) in enumerate(self.elefam[e].edges):
                if tuple(sorted(edge)) == w:
                    # the internal node numbers match, this is the edge
                    surface.append((e, ie))
        return np.array(surface)

    def node_set(self, name, region):
        if not is_stringlike(name):
            raise UserInputError("Name must be a string")
        self.nodesets[name.upper()] = self.get_internal_node_ids(region)

    def side_set(self, name, surface):
        if self.unassigned:
            raise UserInputError(
                "Mesh element operations require all elements be "
                "assigned to an element block"
            )
        if not is_stringlike(name):
            raise UserInputError("Name must be a string")
        self.surfaces[name.upper()] = self.find_surface(surface)

    def element_set(self, name, region):
        if self.unassigned:
            raise UserInputError(
                "Mesh element operations require all elements be "
                "assigned to an element block"
            )
        if not is_stringlike(name):
            raise UserInputError("Name must be a string")
        if region == ALL:
            ielems = range(self.numele)
        else:
            if not is_listlike(region):
                region = [region]
            ielems = [self.elemap[el] for el in region]
        self.elemsets[name.upper()] = np.array(ielems, dtype=int)

    def element_block(self, name, elements):

        if name.upper() in [eb.name.upper() for eb in self.element_blocks]:
            raise UserInputError("{0!r} already an element block".format(name))
        if elements == ALL:
            xelems = sorted(self.elemap.keys())
        else:
            if not is_listlike(elements):
                elements = [elements]
            xelems = elements

        if not self.unassigned:
            raise UserInputError("All elements have been assigned")

        # Elements are numbered in the order the appear in element blocks,
        # we need to map from the temporary internal numbers
        blkcon, badel = [], []
        numblkel = len(xelems)
        ielems = np.arange(self.num_assigned, self.num_assigned + numblkel)
        for (i, xelem) in enumerate(xelems):
            # external element label -> new internal element ID
            ielem = ielems[i]
            self.elemap[xelem] = ielem
            self.ielemap[ielem] = xelem
            elenod = self.eletab[xelem]
            if blkcon and len(elenod) != len(blkcon[-1]):
                badel.append(xelem)
                continue
            blkcon.append(elenod)
        if badel:
            badel = ",\n   ".join(str(e) for e in badel)
            raise UserInputError(
                "The following elements have inconsistent element "
                "connectivity:\n   {0}".format(badel)
            )
        blkcon = np.array(blkcon, dtype=int)
        elefam = element_family(self.dimensions, blkcon.shape[1])
        blk = element_block(name, len(self.element_blocks) + 1, xelems, elefam, blkcon)
        self.element_blocks.append(blk)
        self.num_assigned += len(ielems)
        self.maxedge = max(self.maxedge, len(elefam.edges))
        return blk

    def to_genesis(self, filename):
        exof = exodusii.File(filename, mode="w")
        if not self.element_blocks:
            # put in a single element block
            self.element_block("ElementBlock1", ALL)
        exof.genesis(
            self.nodmap,
            self.elemap,
            self.coord,
            self.element_blocks,
            nodesets=self.nodesets,
            elemsets=self.elemsets,
            sidesets=self.surfaces,
        )
        exof.close()

    def Plot2D(
        self,
        xy=None,
        elecon=None,
        u=None,
        color=None,
        ax=None,
        show=0,
        weight=None,
        colorby=None,
        linestyle="-",
        label=None,
        xlim=None,
        ylim=None,
        filename=None,
        **kwds,
    ):
        assert self.dimensions == 2
        from matplotlib.patches import Polygon
        # import matplotlib.lines as mlines
        from matplotlib.collections import PatchCollection
        from matplotlib.cm import Spectral  # , coolwarm
        import matplotlib.pyplot as plt

        if xy is None:
            xy = np.array(self.coord)
        if elecon is None:
            elecon = []
            for blk in self.element_blocks:
                elecon.extend(blk.elecon.tolist())
            elecon = np.asarray(elecon)
        if u is not None:
            xy += u.reshape(xy.shape)

        patches = []
        for points in xy[elecon[:]]:
            quad = Polygon(points, True)
            patches.append(quad)

        if ax is None:
            fig, ax = plt.subplots()

        # colors = 100 * random.rand(len(patches))
        p = PatchCollection(patches, linewidth=weight, **kwds)
        if colorby is not None:
            colorby = np.asarray(colorby).flatten()
            if len(colorby) == len(xy):
                # average value in element
                colorby = np.array([np.average(colorby[points]) for points in elecon])
            p.set_cmap(Spectral)  # coolwarm)
            p.set_array(colorby)
            p.set_clim(vmin=colorby.min(), vmax=colorby.max())
            fig.colorbar(p)
        else:
            if color is None:
                color = "black"
            p.set_edgecolor(color)
            p.set_facecolor("None")
            p.set_linewidth(weight)
            p.set_linestyle(linestyle)

        if label:
            ax.plot([], [], color=color, linestyle=linestyle, label=label)

        ax.add_collection(p)

        if not ylim:
            ymin, ymax = np.amin(xy[:, 1]), np.amax(xy[:, 1])
            dy = max(abs(ymin * 0.05), abs(ymax * 0.05))
            ax.set_ylim([ymin - dy, ymax + dy])
        else:
            ax.set_ylim(ylim)

        if not xlim:
            xmin, xmax = np.amin(xy[:, 0]), np.amax(xy[:, 0])
            dx = max(abs(xmin * 0.05), abs(xmax * 0.05))
            ax.set_xlim([xmin - dx, xmax + dx])
        else:
            ax.set_xlim(xlim)
        ax.set_aspect("equal")

        if show:
            if label:
                plt.legend()
            plt.show()

        if filename is not None:
            plt.legend()
            plt.savefig(filename, transparent=True, bbox_inches="tight", pad_inches=0)

        return ax

    def PlotScalar2D(self, u, show=0):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        from matplotlib.cm import Spectral

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        elecon = []
        for eb in self.element_blocks:
            elecon.extend(eb.elecon)
        elecon = np.asarray(elecon)
        ax.plot_trisurf(
            self.coord[:, 0], self.coord[:, 1], u, triangles=elecon, cmap=Spectral
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        if show:
            plt.show()
        return

    def put_nodal_solution(self, filename, u):
        import felab.mesh.exodusii as exodusii

        if not self.element_blocks:
            self.element_block("ElementBlock1", ALL)
        if not filename.endswith((".exo", ".e")):
            filename += ".exo"
        exodusii.put_nodal_solution(
            filename, self.nodmap, self.elemap, self.coord, self.element_blocks, u
        )
