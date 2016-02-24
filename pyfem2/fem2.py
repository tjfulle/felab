from numpy import array, dot, zeros
from numpy.linalg import solve, LinAlgError
from fem1 import FiniteElementModel
from elemlibN_2 import ElasticLinknD2
from constants import *

__all__ = ['TrussModel']

class TrussModel(FiniteElementModel):
    numdim = None

    def init(self):
        # Request allocation of field variables
        self.request_output_variable('U', VECTOR, NODE)
        self.request_output_variable('R', VECTOR, NODE)
        self.request_output_variable('S', SCALAR, ELEMENT)
        self.request_output_variable('P', SCALAR, ELEMENT)

    def Solve(self):
        """Solves the finite element problem

        """
        # active DOF set dynamically
        self.active_dof = range(self.elements[0].ndof)
        self.validate(ElasticLinknD2, one=True)
        # Assemble the global stiffness and force
        K = self.assemble_global_stiffness()
        F, Q = self.assemble_global_force()
        Kbc, Fbc = self.apply_bc(K, F+Q)
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')

        # Total force, including reaction, and reaction
        Ft = dot(K, self.dofs)
        R = Ft - F - Q

        # reshape R to be the same shape as coord
        u = self.dofs.reshape(self.numnod, -1)
        R = R.reshape(self.numnod, -1)
        p = self.internal_forces(u)
        s = self.stresses(p)
        self.snapshot(U=u, R=R, P=p, S=s)

    # ----------------------------------------------------------------------- #
    # --------------------------- POSTPROCESSING ---------------------------- #
    # ----------------------------------------------------------------------- #
    def internal_forces(self, u):
        """

        .. _truss_int_force:

        Computes the axial internal forces for all elements in the truss

        Parameters
        ----------
        u : array_like
            nodal displacements
            u[i,j] is the jth coordinate displacement of the ith node

        Returns
        -------
        p : ndarray
            Axial internal force

        See Also
        --------
        pyfem2.elemlibN_2.ElasticLinknD2.internal_force

        """
        p = zeros(self.numele)
        for el in self.elements:
            iel = self.mesh.elemap[el.label]
            p[iel] = el.internal_force(u[el.nodes])
        return p

    def stresses(self, p):
        """

        .. _truss_stresses:

        Computes the stress in each truss element

        Parameters
        ----------
        p : ndarray
            Element axial force
            p[e] is the internal force in the eth element

        Returns
        -------
        s : ndarray
            Array of stresses

        """
        return p / array([el.A for el in self.elements])

# --------------------------------------------------------------------------- #
# -------------------------------- TESTS ------------------------------------ #
# --------------------------------------------------------------------------- #
def test_1():
    from numpy import allclose
    from elemlibN_2 import ElasticLink3D2
    from constants import ALL, X, Y, Z
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
    u = V.steps[-1].field_outputs['U'].data
    R = V.steps[-1].field_outputs['R'].data
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

def test_1a():
    from numpy import allclose
    from elemlibN_2 import ElasticLink2D2
    from constants import ALL, X, Y, Z
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
    u = V.steps[-1].field_outputs['U'].data
    R = V.steps[-1].field_outputs['R'].data
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

def test_2a():
    from constants import ALL, X, Y, Z
    from numpy import allclose
    from elemlibN_2 import ElasticLink2D2
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
    u = V.steps[-1].field_outputs['U'].data
    assert allclose([[-0.05,     0.0882842],
                     [ 0.,       0.,      ],
                     [ 0.,       0.,      ]], u)

def test_2b():
    from constants import ALL, X, Y, Z
    from numpy import allclose
    from elemlibN_2 import ElasticLink3D2
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
    u = V.steps[-1].field_outputs['U'].data
    assert allclose([[-0.05,     0.0882842, 0.],
                     [ 0.,       0.,        0.],
                     [ 0.,       0.,        0.]], u)

def test_3():
    from constants import ALL, X, Y, Z
    from numpy import allclose
    from elemlibN_2 import ElasticLink3D2
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
    u = V.steps[-1].field_outputs['U'].data
    ua = array([[-8.5337228, 0., -31.9486913],
                [ 0.,        0.,   0.       ],
                [ 0.,        0.,   0.       ],
                [ 0.,        0.,   0.       ]])/10
    assert allclose(u, ua)

def test_4():
    from constants import ALL, X, Y, Z
    from numpy import allclose
    from elemlibN_2 import ElasticLink3D2
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
    u = V.steps[-1].field_outputs['U'].data

    assert allclose([[0.00851510679597,0.349956039184,-0.0221277138856],
                     [0.0319156311642,0.349956039184,-0.0322420125936],
                     [0.0115296378344,-0.00976991195147,-0.108526000994],
                     [-0.00403948591251,-0.00878106640481,-0.115392665916],
                     [0.000447704186587,-0.00508916981698,0.0705078974296],
                     [0.00704244773528,-0.00410032427032,0.0773745623519],
                     [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], u)

    V.WriteResults('Job2')

if __name__ == '__main__':
    test_1()
    test_1a()
    test_2a()
    test_2b()
    test_3()
    test_4()
