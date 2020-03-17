import numpy as np

from felab.error import UserInputError
from felab.util.lang import is_stringlike


def plot2d(mesh=None, model=None, **kwds):
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
    if mesh is not None:
        return _plot2d(mesh, **kwds)

    if model is None:
        raise ValueError("model or mesh required")

    deformed = kwds.pop("deformed", False)
    colorby = kwds.pop("colorby", None)
    scale = kwds.pop("scale", 1.0)

    assert model.dimensions == 2
    if model.dimensions != 2:
        raise UserInputError("Plot2D IS ONLY APPLICABLE TO 2D PROBLEMS")
    xy = np.array(model.mesh.coord)
    if deformed:
        xy += scale * model.steps.last.dofs.reshape(xy.shape)
    elecon = []
    for blk in model.mesh.element_blocks:
        if (blk.eletyp.dimensions, blk.eletyp.nodes) == (2, 8):
            raise NotImplementedError("PLOTTING VALID ONLY FOR LINEAR ELEMENT")
        else:
            elecon.extend(blk.elecon)

    if colorby is not None and is_stringlike(colorby):
        colorby = model._get_field(colorby)
    return _plot2d(model.mesh, xy=xy, elecon=np.array(elecon), colorby=colorby, **kwds)


def _plot2d(
    mesh,
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
    assert mesh.dimensions == 2
    from matplotlib.patches import Polygon

    # import matplotlib.lines as mlines
    from matplotlib.collections import PatchCollection
    from matplotlib.cm import Spectral  # , coolwarm
    import matplotlib.pyplot as plt

    if xy is None:
        xy = np.array(mesh.coord)
    if elecon is None:
        elecon = []
        for blk in mesh.element_blocks:
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


def plot2d_scalar(model, u, show=0):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from matplotlib.cm import Spectral

    mesh = model.mesh

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    elecon = []
    for eb in mesh.element_blocks:
        elecon.extend(eb.elecon)
    elecon = np.asarray(elecon)
    ax.plot_trisurf(
        mesh.coord[:, 0], mesh.coord[:, 1], u, triangles=elecon, cmap=Spectral
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if show:
        plt.show()
    return
