from numpy import *
from numpy.linalg import solve, LinAlgError
from constants import *
from utilities import *
from fem1 import FiniteElementModel
from elemlib2_3T import DiffussiveHeatTransfer2D3

class HeatTransfer2DModel(FiniteElementModel):
    numdim = 2
    def init(self):
        # Request allocation of field variables
        self.request_output_variable('T', SCALAR, NODE)
        self.request_output_variable('R', SCALAR, NODE)

    def Solve(self):
        self.validate(DiffussiveHeatTransfer2D3)
        K = self.assemble_global_stiffness(self.sfilm)
        F, Q = self.assemble_global_force(self.src, self.sflux, self.sfilm)
        Kbc, Fbc = self.apply_bc(K, F+Q)
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')
        Ft = dot(K, self.dofs)
        R = Ft - F - Q
        self.snapshot(T=self.dofs, R=R)
        self.T = self.dofs

def plate_with_hole():
    from elemlib2_3T import DiffussiveHeatTransfer2D3
    V = HeatTransfer2DModel()
    V.GenesisMesh('../meshes/PlateWithHoleTria3Fine.g')
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
    V.mesh.PlotScalar2D(V.dofs.flatten(), show=1)

def test_1():
    from numpy import sin, sinh, mean
    from elemlib2_3T import DiffussiveHeatTransfer2D3
    def solution(x, N=20):
        def fun(k):
            a = sin(k * pi * (1 + x[:,0]) / 2.) / (k ** 3 * sinh(k * pi))
            b = sinh(k * pi * (1 + x[:,1]) / 2.) + sinh(k * pi * (1 - x[:,1]) / 2.)
            return a * b
        u = (1-x[:,0]**2)/2.-16./pi**3*sum([fun(k) for k in range(1, N, 2)],0)
        return u
    V = HeatTransfer2DModel()
    V.GenesisMesh('../meshes/UniformPlateTria3Fine.g')
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.HeatGeneration(ALL, 1)
    V.PrescribedBC(BOUNDARY, T, 0)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 1e-4

def test_2():
    from numpy import mean
    from elemlib2_3T import DiffussiveHeatTransfer2D3
    def solution(x):
        return 2. * (1. + x[:,1]) / ((3. + x[:,0])**2 + (1 + x[:,1])**2)
    V = HeatTransfer2DModel()
    V.GenesisMesh('../meshes/UniformPlateTria3.g')
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(1.)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(BOUNDARY, T, solution)
    V.HeatGeneration(ALL, 0)
    V.Solve()
    Tn = solution(V.mesh.coord)
    err = sqrt(mean((V.dofs.flatten()-Tn)**2)) / sqrt(mean(Tn**2))
    assert err < 5e-5

def test_3():
    from numpy import mean, random, where, sin, sinh
    from distmesh import drectangle, distmesh2d, huniform
    from elemlib2_3T import DiffussiveHeatTransfer2D3
    def solution(x, N=20):
        def fun(n):
            return sin(2.*n*pi/3.)/n**2/sinh(n*pi)*sin(n*pi*x[:,0])*sinh(n*pi*x[:,1])
        return 450./pi**2*sum([fun(n) for n in range(1, N)], 0)
    random.seed(190) # Always the same results
    fd = lambda p: drectangle(p, 0, 1, 0, 1)
    fh = huniform
    coord, elecon = distmesh2d(fd, fh, .05, (0, 0, 1, 1),
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

def test_4():
    from numpy import cos, cosh, mean
    from elemlib2_3T import DiffussiveHeatTransfer2D3
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

def test_5():
    from numpy import mean, random, where, sin, sinh
    from distmesh import drectangle, distmesh2d, huniform
    from numpy import cos, cosh, mean
    from elemlib2_3T import DiffussiveHeatTransfer2D3
    def solution(x, q0=1., k=1., N=20):
        def fun(n):
            al = .5 * (2. * n - 1.) * pi
            top = (-1)**n*cos(al*x[:,1])*cosh(al*x[:,0])
            bot = al**3*cosh(al)
            return top / bot
        return q0/2/k*((1-x[:,1]**2)+4.*sum([fun(n) for n in range(1,N)], 0))
    random.seed(190) # Always the same results
    fd = lambda p: drectangle(p, 0, 1, 0, 1)
    fh = huniform
    coord, elecon = distmesh2d(fd, fh, .025, (0, 0, 1, 1),
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

def plate_with_hole2():
    from elemlib2_3T import DiffussiveHeatTransfer2D3
    k, h, Too = 12, 250, 25
    V = HeatTransfer2DModel()
    V.GenesisMesh('../meshes/PlateWithHoleTria3.g')
    V.Material('Material-1')
    V.materials['Material-1'].IsotropicThermalConductivity(k)
    V.AssignProperties('ElementBlock1', DiffussiveHeatTransfer2D3, 'Material-1')
    V.PrescribedBC(BOUNDARY, T, 50)
    fun = lambda x: 1000 / sqrt(x[:,0] ** 2 + x[:,1] ** 2)
    V.HeatGeneration(ALL, fun)
    V.InitialTemperature(ALL, 50)
    V.Solve()
    V.WriteResults('Job3')

if __name__ == '__main__':
    test_1()
    test_2()
    test_3()
    test_4()
    test_5()
    plate_with_hole2()
