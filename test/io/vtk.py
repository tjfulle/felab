import pytest
try:
    import distmesh
except ImportError:
    distmesh = None
import felab.util.tty as tty
from felab.io.vtk import WriteFEResults, WriteVTUMesh


@pytest.mark.skipif(distmesh is None, reason="Requires distmesh")
def test_write_fe_results():
    if distmesh is None:
        tty.warn("test requires distmesh module")
        return
    np.random.seed(190)  # Always the same results
    fd = lambda p: dm.drectangle(p, -1, 1, -1, 1)
    fh = dm.huniform
    coord, elecon = dm.distmesh2d(
        fd, fh, 0.1, (-1, -1, 1, 1), [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    )
    jobid = "Job"
    nodlab = range(coord.shape[0])
    nodmap = dict([(n, n) for n in nodlab])
    elelab = range(elecon.shape[0])
    elemap = dict([(n, n) for n in elelab])
    eletyp = [element_family(2, 3)] * elecon.shape[0]
    scal = np.random.rand(coord.shape[0])
    vect = np.random.rand(coord.shape[0] * 2).reshape(-1, 2)
    tens = np.random.rand(elecon.shape[0] * 9).reshape(-1, 9)
    symt = np.random.rand(elecon.shape[0] * 6).reshape(-1, 6)
    kwds = dict(scal=scal, vect=vect, tens=tens, symt=symt)
    u = np.zeros_like(coord)
    u[:, 0] = 1
    WriteFEResults(jobid, coord, nodmap, elemap, eletyp, elecon, u=u, **kwds)
    filename = jobid + ".vtu"
    WriteVTUMesh(filename, coord, nodlab, elelab, eletyp, elecon, check=1)
    os.remove(filename)
