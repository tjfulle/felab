from numpy import *
from conf import *
from felab import *


def test_gravity_load1():
    V = fe_model(jobid='Gravity')
    V.rectilinear_mesh(nx=101, ny=11, lx=100, ly=10)
    V.create_material('Material-1')
    V.materials['Material-1'].Density(1.)
    V.materials['Material-1'].elastic(E=10e6, Nu=.333)
    V.create_element_block('Block1', ALL)
    V.assign_properties('Block1', CPE4, 'Material-1', t=1)
    stage = V.create_static_stage()
    stage.fix_nodes(ILO)
    stage.assign_gravity_load(ALL, [0, -1e2])
    stage.run()
    V.write_results()

