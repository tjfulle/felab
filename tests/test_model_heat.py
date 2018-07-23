import pytest
from numpy import *
from conf import *
try:
    import distmesh as dm
except ImportError:
    dm = None
from felab import *

def test_heat_transfer_1():
    def solution(x, N=20):
        def fun(k):
            a = sin(k * pi * (1 + x[:,0]) / 2.) / (k ** 3 * sinh(k * pi))
            b = sinh(k * pi * (1 + x[:,1]) / 2.) + sinh(k * pi * (1 - x[:,1]) / 2.)
            return a * b
        u = (1-x[:,0]**2)/2.-16./pi**3*sum([fun(k) for k in range(1, N, 2)],0)
        return u
    V = fe_model()
    V.genesis_mesh('./data/UniformPlateTria3Fine.g')
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(1.)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    stage = V.create_heat_transfer_stage()
    stage.HeatGeneration(ALL, 1)
    stage.assign_prescribed_bc(BOUNDARY, T, 0)
    stage.run()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((stage.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 1e-4

def test_heat_transfer_2():
    def solution(x):
        return 2. * (1. + x[:,1]) / ((3. + x[:,0])**2 + (1 + x[:,1])**2)
    V = fe_model()
    V.genesis_mesh('./data/UniformPlateTria3.g')
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(1.)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(BOUNDARY, T, solution)
    stage.HeatGeneration(ALL, 0)
    stage.run()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((stage.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 5e-5

@pytest.mark.skipif(dm is None, reason='distmesh not import')
def test_heat_transfer_3():
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
    V = fe_model()
    V.create_mesh(p=coord, t=elecon)
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(1.)
    V.create_element_block('ElementBlock1', ALL)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(JLO, T, 0)
    stage.assign_prescribed_bc(JHI, T, f2)
    stage.assign_prescribed_bc(ILO, T, 0)
    stage.assign_prescribed_bc(IHI, T, 0)
    stage.HeatGeneration(ALL, 0)
    stage.run()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((stage.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 5e-3

def test_heat_transfer_3():
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
    V = fe_model()
    V.create_mesh(nodtab=nodtab, eletab=eletab)
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(1.)
    V.create_element_block('ElementBlock1', ALL)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(IHI, T)
    stage.assign_prescribed_bc(JHI, T)
    stage.HeatGeneration(ALL, 1)
    stage.run()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((stage.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < .045

@pytest.mark.skipif(dm is None, reason='distmesh not import')
def test_heat_transfer_4():
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
    V = fe_model()
    V.create_mesh(p=coord, t=elecon)
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(1.)
    V.create_element_block('ElementBlock1', ALL)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(IHI, T)
    stage.assign_prescribed_bc(JHI, T)
    stage.HeatGeneration(ALL, 1)
    stage.run()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((stage.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 1e-4

def test_heat_transfer_plate_with_hole_1():
    V = fe_model()
    V.genesis_mesh('./data/PlateWithHoleTria3Fine.g')
    k, h, Too = 12, 250, 25
    fun = lambda x: 1000 / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(k)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(ILO, T, 200)
    stage.assign_prescribed_bc(IHI, T, 50)
    stage.SurfaceFlux(JLO, 2000)
    stage.SurfaceConvection(JHI, Too, h)
    stage.HeatGeneration(ALL, fun)
    stage.run()
    V.mesh.PlotScalar2D(stage.dofs.flatten())

def test_heat_transfer_plate_with_hole2():
    k, h, Too = 12, 250, 25
    V = fe_model(jobid='HeatPlateWithHole')
    V.genesis_mesh('./data/PlateWithHoleTria3.g')
    V.create_material('Material-1')
    V.materials['Material-1'].isotropic_thermal_conductivity(k)
    V.assign_properties('ElementBlock1', DC2D3, 'Material-1')
    V.assign_initial_temperature(ALL, 50)
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(BOUNDARY, T, 50)
    fun = lambda x: 1000 / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
    stage.HeatGeneration(ALL, fun)
    stage.run()
    V.write_results()

