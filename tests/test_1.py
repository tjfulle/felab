import os
import re
import sys
import glob
import pytest
import shutil
from subprocess import Popen
from os.path import join, dirname, realpath, isfile, isdir, basename
from numpy import cos, cosh, mean, random, where, sin, sinh, allclose
from conf import *
try:
    import distmesh as dm
except ImportError:
    dm = None

PY3 = sys.version_info[0] == 3
D = os.getenv('PYFEM2')
if D is None:
    D = dirname(dirname(realpath(__file__)))
sys.path.insert(0, D)
from pyfem2 import *

# --------------------------------------------------------------------------- #
# --------------------------------- TESTS ----------------------------------- #
# --------------------------------------------------------------------------- #
@pytest.mark.beamcol
def test_fem5_1():
    nodtab = [[1,-4,3], [2,0,0], [3,0,3], [4,nan,nan], [5,4,3]]
    eletab = [[1,1,3], [2,3,5], [3,1,2], [4,2,3], [5,2,5]]
    V = PlaneBeamColumnTrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    Ec, Em = 30000, 200000
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=Ec, Nu=.3)
    V.Material('Material-2')
    V.materials['Material-2'].Elastic(E=Em, Nu=.3)
    V.ElementBlock('B1', (1,2))
    V.ElementBlock('B2', (3,5))
    V.ElementBlock('B3', (4,))
    V.AssignProperties('B1', BeamColumn2D, 'Material-1', A=.02, Izz=.004)
    V.AssignProperties('B2', ElasticLink2D2, 'Material-2', A=.001)
    V.AssignProperties('B3', ElasticLink2D2, 'Material-2', A=.003)
    V.PrescribedBC(1, (X,Y,TZ))
    V.PrescribedBC(5, Y)
    V.ConcentratedLoad(2, Y, 100)
    V.ConcentratedLoad(5, TZ, 200)
    V.ConcentratedLoad(5, X, 300)
    V.Solve()

@pytest.mark.beamcol
def test_fem5_2():
    nodtab = [[1,0,0], [2,3,4], [3,0,4]]
    eletab = [[1,1,2], [2,1,3]]
    V = PlaneBeamColumnTrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=70e9, Nu=.3)
    V.ElementBlock('B1', ALL)
    V.AssignProperties('B1', ElasticLink2D2, 'Material-1', A=5*.01*.01)
    V.PrescribedBC(1, X, -.05)
    V.PrescribedBC((2,3), (X,Y))
    V.ConcentratedLoad(1, Y, 1000e3)
    V.Solve()
    assert allclose([[-0.05,     0.0882842],
                     [ 0.,       0.,      ],
                     [ 0.,       0.,      ]], V.u)

@pytest.mark.plane
def test_fem4_1_gravity_load1():
    V = Plane2DModel('Gravity')
    V.RectilinearMesh((101, 11), (100, 10))
    V.Material('Material-1')
    V.materials['Material-1'].Density(1.)
    V.materials['Material-1'].Elastic(E=10e6, Nu=.333)
    V.ElementBlock('Block1', ALL)
    V.AssignProperties('Block1', PlaneStrainQuad4, 'Material-1', t=1)
    V.FixNodes(ILO)
    V.GravityLoad(ALL, [0, -1e2])
    V.Solve()
    V.WriteResults()

@pytest.mark.heat
def test_fem3_plate_with_hole():
    V = HeatTransfer2DModel()
    V.GenesisMesh(join(D, 'data/PlateWithHoleTria3Fine.g'))
    k, h, Too = 12, 250, 25
    fun = lambda x: 1000 / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(k)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(ILO, T, 200)
    V.PrescribedBC(IHI, T, 50)
    V.SurfaceFlux(JLO, 2000)
    V.SurfaceConvection(JHI, Too, h)
    V.HeatGeneration(ALL, fun)
    V.Solve()
    V.mesh.PlotScalar2D(V.dofs.flatten())

@pytest.mark.heat
def test_fem3_1():
    def solution(x, N=20):
        def fun(k):
            a = sin(k * pi * (1 + x[:,0]) / 2.) / (k ** 3 * sinh(k * pi))
            b = sinh(k * pi * (1 + x[:,1]) / 2.) + sinh(k * pi * (1 - x[:,1]) / 2.)
            return a * b
        u = (1-x[:,0]**2)/2.-16./pi**3*sum([fun(k) for k in range(1, N, 2)],0)
        return u
    V = HeatTransfer2DModel()
    V.GenesisMesh(join(D, 'data/UniformPlateTria3Fine.g'))
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.HeatGeneration(ALL, 1)
    V.PrescribedBC(BOUNDARY, T, 0)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 1e-4

@pytest.mark.heat
def test_fem3_2():
    def solution(x):
        return 2. * (1. + x[:,1]) / ((3. + x[:,0])**2 + (1 + x[:,1])**2)
    V = HeatTransfer2DModel()
    V.GenesisMesh(join(D, 'data/UniformPlateTria3.g'))
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(BOUNDARY, T, solution)
    V.HeatGeneration(ALL, 0)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 5e-5

@pytest.mark.heat
@pytest.mark.skipif(dm is None, reason='distmesh not import')
def test_fem3_3():
    def solution(x, N=20):
        def fun(n):
            return sin(2.*n*pi/3.)/n**2/sinh(n*pi)*sin(n*pi*x[:,0])*sinh(n*pi*x[:,1])
        return 450./pi**2*sum([fun(n) for n in range(1, N)], 0)
    random.seed(190) # Always the same results
    fd = lambda p: dm.drectangle(p, 0, 1, 0, 1)
    fh = dm.huniform
    coord, elecon = dm.distmesh2d(fd, fh, .05, (0, 0, 1, 1),
                                  [(0, 0),(0, 1),(1, 0),(1, 1)])
    f2 = lambda x: where(x[:,0] <= 2./3., 75*x[:,0], 150*(1-x[:,0]))
    V = HeatTransfer2DModel()
    V.Mesh(p=coord, t=elecon)
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(JLO, T, 0)
    V.PrescribedBC(JHI, T, f2)
    V.PrescribedBC(ILO, T, 0)
    V.PrescribedBC(IHI, T, 0)
    V.HeatGeneration(ALL, 0)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 5e-3

@pytest.mark.heat
def test_fem3_4():
    def solution(x, q0=1., k=1., N=20):
        def fun(n):
            al = .5 * (2. * n - 1.) * pi
            top = (-1)**n*cos(al*x[:,1])*cosh(al*x[:,0])
            bot = al**3*cosh(al)
            return top / bot
        return q0/2/k*((1-x[:,1]**2)+4.*sum([fun(n) for n in range(1,N)], 0))
    nodtab = [[1,0.,0.], [2,.5,0.], [3,1.,0.],
              [4,0.,.5], [5,.5,.5], [6,1.,.5],
              [7,0.,1.], [8,.5,1.], [9,1.,1.]]
    eletab = [[1,1,2,5], [2,1,5,4], [3,2,3,6], [4,2,6,5],
              [5,4,5,8], [6,4,8,7], [7,5,6,9], [8,5,9,8]]
    V = HeatTransfer2DModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(IHI, T)
    V.PrescribedBC(JHI, T)
    V.HeatGeneration(ALL, 1)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < .045

@pytest.mark.heat
@pytest.mark.skipif(dm is None, reason='distmesh not import')
def test_fem3_5():
    def solution(x, q0=1., k=1., N=20):
        def fun(n):
            al = .5 * (2. * n - 1.) * pi
            top = (-1)**n*cos(al*x[:,1])*cosh(al*x[:,0])
            bot = al**3*cosh(al)
            return top / bot
        return q0/2/k*((1-x[:,1]**2)+4.*sum([fun(n) for n in range(1,N)], 0))
    random.seed(190) # Always the same results
    fd = lambda p: dm.drectangle(p, 0, 1, 0, 1)
    fh = dm.huniform
    coord, elecon = dm.distmesh2d(fd, fh, .025, (0, 0, 1, 1),
                                  [(0, 0),(0, 1),(1, 0),(1, 1)])
    V = HeatTransfer2DModel()
    V.Mesh(p=coord, t=elecon)
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(IHI, T)
    V.PrescribedBC(JHI, T)
    V.HeatGeneration(ALL, 1)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 1e-4

@pytest.mark.heat
def test_fem3_plate_with_hole2():
    k, h, Too = 12, 250, 25
    V = HeatTransfer2DModel('HeatPlateWithHole')
    V.GenesisMesh(join(D, 'data/PlateWithHoleTria3.g'))
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(k)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(BOUNDARY, T, 50)
    fun = lambda x: 1000 / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
    V.HeatGeneration(ALL, fun)
    V.InitialTemperature(ALL, 50)
    V.Solve()
    V.WriteResults()

@pytest.mark.truss
def test_fem2_1():
    nodtab = [[1,0,0,0], [2,10,5,0], [3,10,0,0], [4,20,8,0], [5,20,0,0],
              [6,30,9,0], [7,30,0,0], [8,40,8,0], [9,40,0,0], [10,50,5,0],
              [11,50,0,0], [12,60,0,0]]
    eletab = [[1,1,3], [2,3,5], [3,5,7], [4,7,9], [5,9,11], [6,11,12],
              [7,1,2], [8,2,4], [9,4,6], [10,6,8], [11,8,10], [12,10,12],
              [13,2,3], [14,4,5], [15,6,7], [16,8,9], [17,10,11], [18,2,5],
              [19,4,7], [20,7,8], [21,9,10]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1000, Nu=.333)
    Abot, Atop, Abat, Adia = 2, 10, 3, 1
    A = [Abot, Abot, Abot, Abot, Abot, Abot,
         Atop, Atop, Atop, Atop, Atop, Atop,
         Abat, Abat, Abat, Abat, Abat,
         Adia, Adia, Adia, Adia]
    V.ElementBlock('ElementBlock-1', ALL)
    V.AssignProperties('ElementBlock-1', ElasticLink3D2, 'Material-1', A=A)

    V.PrescribedBC(1, (X,Y))
    V.PrescribedBC(12, Y)
    V.PrescribedBC(ALL, Z)

    V.ConcentratedLoad((3,5,9,11), Y, -10)
    V.ConcentratedLoad(7, Y, -16)

    V.Solve()
    u = V.steps.last.frames[-1].field_outputs['U'].data
    R = V.steps.last.frames[-1].field_outputs['R'].data
    assert allclose([[0.,      0.,      0.],
                     [0.80954,-1.7756,  0.],
                     [0.28,   -1.79226, 0.],
                     [0.899,  -2.29193, 0.],
                     [0.56,   -2.3166,  0.],
                     [0.8475, -2.38594, 0.],
                     [0.8475, -2.42194, 0.],
                     [0.796,  -2.29193, 0.],
                     [1.135,  -2.3166,  0.],
                     [0.88546,-1.7756,  0.],
                     [1.415,  -1.79226, 0.],
                     [1.695,   0.,      0.]], u)
    assert allclose([[ 0.,  28.,   0.],
                     [ 0.,   0.,   0.],
                     [ 0.,   0.,   0.],
                     [ 0.,   0.,   0.],
                     [ 0.,   0.,   0.],
                     [-0.,  -0.,   0.],
                     [-0.,   0.,   0.],
                     [-0.,   0.,   0.],
                     [ 0.,  -0.,   0.],
                     [ 0.,   0.,   0.],
                     [ 0.,   0.,   0.],
                     [ 0.,  28.,   0.]], R)

@pytest.mark.truss
def test_fem2_1a():
    nodtab = [[1,0,0], [2,10,5], [3,10,0], [4,20,8], [5,20,0],
              [6,30,9], [7,30,0], [8,40,8], [9,40,0], [10,50,5],
              [11,50,0], [12,60,0]]
    eletab = [[1,1,3], [2,3,5], [3,5,7], [4,7,9], [5,9,11], [6,11,12],
              [7,1,2], [8,2,4], [9,4,6], [10,6,8], [11,8,10], [12,10,12],
              [13,2,3], [14,4,5], [15,6,7], [16,8,9], [17,10,11], [18,2,5],
              [19,4,7], [20,7,8], [21,9,10]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=1000, Nu=.333)
    Abot, Atop, Abat, Adia = 2, 10, 3, 1
    A = [Abot, Abot, Abot, Abot, Abot, Abot,
         Atop, Atop, Atop, Atop, Atop, Atop,
         Abat, Abat, Abat, Abat, Abat,
         Adia, Adia, Adia, Adia]
    V.ElementBlock('ElementBlock-1', ALL)
    V.AssignProperties('ElementBlock-1', ElasticLink2D2, 'Material-1', A=A)

    V.PrescribedBC(1, (X,Y))
    V.PrescribedBC(12, Y)

    V.ConcentratedLoad((3,5,9,11), Y, -10)
    V.ConcentratedLoad(7, Y, -16)

    V.Solve()
    u = V.steps.last.frames[-1].field_outputs['U'].data
    R = V.steps.last.frames[-1].field_outputs['R'].data
    assert allclose([[0.,      0.     ],
                     [0.80954,-1.7756 ],
                     [0.28,   -1.79226],
                     [0.899,  -2.29193],
                     [0.56,   -2.3166 ],
                     [0.8475, -2.38594],
                     [0.8475, -2.42194],
                     [0.796,  -2.29193],
                     [1.135,  -2.3166 ],
                     [0.88546,-1.7756 ],
                     [1.415,  -1.79226],
                     [1.695,   0.     ]], u)
    assert allclose([[ 0.,  28.],
                     [ 0.,   0.],
                     [ 0.,   0.],
                     [ 0.,   0.],
                     [ 0.,   0.],
                     [-0.,  -0.],
                     [-0.,   0.],
                     [-0.,   0.],
                     [ 0.,  -0.],
                     [ 0.,   0.],
                     [ 0.,   0.],
                     [ 0.,  28.]], R)

@pytest.mark.truss
def test_fem2_2a():
    nodtab = [[1, 0., 0.], [2, 3., 4.], [3, 0., 4.]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=70e9, Nu=.333)
    A = 5 * .01 * .01
    V.ElementBlock('ElementBlock-1', ALL)
    V.AssignProperties('ElementBlock-1', ElasticLink2D2, 'Material-1', A=A)
    V.FixNodes((2,3))
    V.PrescribedBC(1, X, -.05)
    V.ConcentratedLoad(1, Y, 1000e3)
    V.Solve()
    u = V.steps.last.frames[-1].field_outputs['U'].data
    assert allclose([[-0.05,     0.0882842],
                     [ 0.,       0.,      ],
                     [ 0.,       0.,      ]], u)

@pytest.mark.truss
def test_fem2_2b():
    nodtab = [[1, 0., 0., 0.], [2, 3., 4., 0.], [3, 0., 4., 0.]]
    eletab = [[1, 1, 2], [2, 1, 3]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=70e9, Nu=.333)
    A = 5 * .01 * .01
    V.ElementBlock('ElementBlock-1', ALL)
    V.AssignProperties('ElementBlock-1', ElasticLink3D2, 'Material-1', A=A)
    V.FixNodes((2,3))
    V.PrescribedBC(1, X, -.05)
    V.PrescribedBC(1, Z)
    V.ConcentratedLoad(1, Y, 1000e3)
    V.Solve()
    u = V.steps.last.frames[-1].field_outputs['U'].data
    assert allclose([[-0.05,     0.0882842, 0.],
                     [ 0.,       0.,        0.],
                     [ 0.,       0.,        0.]], u)

@pytest.mark.truss
def test_fem2_3():
    nodtab = [[1, 72, 0, 0], [2, 0, 36, 0], [3, 0, 36, 72], [4, 0, 0, -48]]
    eletab = [[1, 1, 2], [2, 1, 3], [3, 1, 4]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=10e4, Nu=.333)
    A = [.302, .729, .187]
    V.ElementBlock('ElementBlock-1', ALL)
    V.AssignProperties('ElementBlock-1', ElasticLink3D2, 'Material-1', A=A)
    # Boundary conditions
    V.FixNodes((2,3,4))
    V.PrescribedBC(1, Y)
    # Concentrated force in 'z' direction on node 1
    V.ConcentratedLoad(1, Z, -1000)
    V.Solve()
    u = V.steps.last.frames[-1].field_outputs['U'].data
    ua = array([[-8.5337228, 0., -31.9486913],
                [ 0.,        0.,   0.       ],
                [ 0.,        0.,   0.       ],
                [ 0.,        0.,   0.       ]])/10
    assert allclose(u, ua)

@pytest.mark.truss
def test_fem2_4():
    # Set up problem space
    nodtab = [[1,-37.5,0,200],[2,37.5,0,200],[3,-37.5,37.5,100],
              [4,37.5,37.5,100],[5,37.5,-37.5,100],[6,-37.5,-37.5,100],
              [7,-100,100,0],[8,100,100,0],[9,100,-100,0],[10,-100,-100,0]]
    eletab = [[ 1, 1, 2],[ 2, 1, 4],[ 3, 2, 3],[ 4, 1, 5],[ 5, 2, 6],
              [ 6, 2, 4],[ 7, 2, 5],[ 8, 1, 3],[ 9, 1, 6],[10, 3, 6],
              [11, 4, 5],[12, 3, 4],[13, 5, 6],[14, 3,10],[15, 6, 7],
              [16, 4, 9],[17, 5, 8],[18, 4, 7],[19, 3, 8],[20, 5,10],
              [21, 6, 9],[22, 6,10],[23, 3, 7],[24, 5, 9],[25, 4, 8]]
    V = TrussModel()
    V.Mesh(nodtab=nodtab, eletab=eletab)

    # Define element blocks
    V.Material('Material-1')
    V.materials['Material-1'].Elastic(E=10e6, Nu=.333)
    A = [0.033, 2.015, 2.015, 2.015, 2.015, 2.823, 2.823, 2.823, 2.823, 0.01,
         0.01, 0.014, 0.014, 0.98, 0.98, 0.98, 0.98, 1.76, 1.76, 1.76, 1.76,
         2.44, 2.44, 2.44, 2.44]
    V.ElementBlock('ElementBlock-1', ALL)
    V.AssignProperties('ElementBlock-1', ElasticLink3D2, 'Material-1', A=A)

    # Define boundary conditons
    V.FixNodes([7, 8, 9, 10])

    # Define concentrated loads
    V.ConcentratedLoad(1, X, 1000)
    V.ConcentratedLoad(1, Y, 10000)
    V.ConcentratedLoad(1, Z, -5000)
    V.ConcentratedLoad(2, Y, 10000)
    V.ConcentratedLoad(2, Z, -5000)
    V.ConcentratedLoad(3, X, 500)
    V.ConcentratedLoad(6, X, 500)

    # Solve and write results
    V.Solve()
    u = V.steps.last.frames[-1].field_outputs['U'].data

    assert allclose([[0.00851510679597,0.349956039184,-0.0221277138856],
                     [0.0319156311642,0.349956039184,-0.0322420125936],
                     [0.0115296378344,-0.00976991195147,-0.108526000994],
                     [-0.00403948591251,-0.00878106640481,-0.115392665916],
                     [0.000447704186587,-0.00508916981698,0.0705078974296],
                     [0.00704244773528,-0.00410032427032,0.0773745623519],
                     [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], u)

def test_fem1_1():
    xa, xb = 0., 1.
    A, E, p, q = 1, 1, 1, 0
    u, e, f = UniformBar(xa, xb, A, E, p, q, numele=10)
    assert abs(u[-1] - 1.) < 1e-12

def test_fem1_2():
    xa, xb = 0., 1.
    A, E, p, q = 1, 1, 1, 1
    u, e, f = UniformBar(xa, xb, A, E, p, q, numele=1000)
    assert abs(u[-1] - 1.5) < 5e-4

def test_elemlibN_2_1():
    """ElasticLink2 stiffness test"""
    class mat: E = 1
    El = ElasticLink1D2(1, [0, 1], [0, 1], mat, A=1)
    Ke = El.response(zeros(2), zeros(2), [0,0], 1., 1, 1, [], [], [], STATIC,
                     [], False, STIFF_ONLY, LINEAR_PERTURBATION, 1.)
    assert allclose([[1,-1],[-1,1]], Ke)

def test_elemlibN_2_2():
    """ElasticLink2 stiffness test"""
    class mat: E=1000
    El = ElasticLink2D2(1, [0, 1], [[0,0], [30,40]], mat, A=5)
    Ke = El.response(zeros((2,2)),zeros((2,2)),[0,0],1.,1,1,[],[],[],STATIC,
                     [], False, STIFF_ONLY, LINEAR_PERTURBATION, 1.)
    assert allclose([[ 36.,  48., -36., -48.],
                     [ 48.,  64., -48., -64.],
                     [-36., -48.,  36.,  48.],
                     [-48., -64.,  48.,  64.]], Ke)

def test_elemlibN_2_3():
    """ElasticLink2 stiffness test"""
    class mat: E = 343
    El = ElasticLink3D2(1, [0, 1], [[0,0,0],[2,3,6]], mat, A=10)
    Ke = El.response(zeros((2,3)),zeros((2,3)),[0,0],1.,1,1,[],[],[],STATIC,
                     [], False, STIFF_ONLY, LINEAR_PERTURBATION, 1.)
    assert allclose([[  40.,   60.,  120.,  -40.,  -60., -120.],
                     [  60.,   90.,  180.,  -60.,  -90., -180.],
                     [ 120.,  180.,  360., -120., -180., -360.],
                     [ -40.,  -60., -120.,   40.,   60.,  120.],
                     [ -60.,  -90., -180.,   60.,   90.,  180.],
                     [-120., -180., -360.,  120.,  180.,  360.]], Ke)

def test_elemlibN_2_4():
    """Beam-Column stiffness test"""
    coord = array([[0, 0], [3, 4]], dtype=float)
    class mat: E = 100
    A, Izz = 125, 250
    El = BeamColumn2D(1, [0, 1], coord, mat, A=A, Izz=Izz)
    Ke = El.response(zeros((2,2)),zeros((2,2)),[0,0],1.,1,1,[],[],[],STATIC,
                     [], False, STIFF_ONLY, LINEAR_PERTURBATION, 1.)

@pytest.mark.demos
def test_demos():
    d = realpath(join(D, 'data'))
    env = dict(os.environ)
    env['PYTHONPATH'] = D
    env['NOGRAPHICS'] = '1'
    failed = []
    cwd = os.getcwd()
    os.chdir(d)
    for filename in glob.glob(join(d, '*.py')):
        if not re.search(r'(?i)Demo[a-z0-9]+\.py', filename):
            continue
        proc = Popen(['python', filename], env=env)
        proc.wait()
        if proc.returncode != 0:
            failed.append(filename)
    os.chdir(cwd)
    assert not failed

if __name__ == '__main__':
    test_demos()
