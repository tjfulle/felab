""".pft.py: Plane Beam Column Truss application code.

"""
from numpy import *
from numpy.linalg import solve, LinAlgError

from constants import *
from fem1 import FiniteElementModel
from elemlibN_2 import ElasticLink2D2, BeamColumn2D

# --------------------------------------------------------------------------- #
# -------------------------- APPLICATION CODE ------------------------------- #
# --------------------------------------------------------------------------- #
class PlaneBeamColumnTrussModel(FiniteElementModel):
    numdim = 2
    def init(self):
        # Request allocation of field variables
        self.request_output_variable('U', VECTOR, NODE)
        self.request_output_variable('Rot', SCALAR, NODE)
        self.request_output_variable('R', VECTOR, NODE)
        self.request_output_variable('Rt', SCALAR, NODE)

    def Solve(self):
        # active DOF set dynamically
        self.validate((ElasticLink2D2, BeamColumn2D))

        # Assemble the global stiffness and force
        K = self.assemble_global_stiffness()
        F, Q = self.assemble_global_force()
        Kbc, Fbc = self.apply_bc(K, F+Q)

        # Solve
        try:
            self.dofs = solve(Kbc, Fbc)
        except LinAlgError:
            raise RuntimeError('attempting to solve under constrained system')

        # Total force, including reaction, and reaction
        Ft = dot(K, self.dofs)
        R = Ft - F - Q

        U, R, Urot, Rrot = self.format_displacements_and_reactions(self.dofs, R)

        self.u = U

        self.snapshot(U=U, R=R, Rot=Urot, Rt=Rrot)

    def format_displacements_and_reactions(self, u, R):
        # construct displacement and rotation vectors
        ut, rt = zeros((self.numnod, 2)), zeros((self.numnod, 2))
        ur, rr = zeros(self.numnod), zeros(self.numnod)
        for n in range(self.numnod):
            ix = 0
            for j in range(MDOF):
                if self.nodfat[n,j] > 0:
                    ii = self.nodfmt[n] + ix
                    if j <= 2:
                        ut[n,j] = u[ii]
                        rt[n,j] = R[ii]
                    else:
                        # This would need modifying if 3D beams were allowed
                        ur[n] = u[ii]
                        rr[n] = R[ii]
                    ix += 1
        return ut, rt, ur, rr

# --------------------------------------------------------------------------- #
# --------------------------------- TESTS ----------------------------------- #
# --------------------------------------------------------------------------- #
def test_1():
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

def test_2():
    from mesh import Mesh
    from mat import Material
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

if __name__ == '__main__':
    test_1()
    test_2()
