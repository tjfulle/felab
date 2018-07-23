from felab import *

def QuarterPlateDemo(plot=False):
    V = fe_model(jobid='PlateWithHoleQuad4QuarterSym')
    V.genesis_mesh('./data/PlateWithHoleQuad4QuarterSym.g')
    mat = V.create_material('Material-1')
    mat.elastic(E=100, Nu=.2)
    V.assign_properties('', CPE4, mat, t=1)
    V.assign_prescribed_bc('SymYZ', X)
    V.assign_prescribed_bc('SymXZ', Y)
    V.assign_initial_temperature(ALL, 60)

    stage = V.create_static_stage('Stage-1')
    stage.assign_surface_load('RightHandSide', [1,0])
    stage.run()

    V.write_results()

    F = File('PlateWithHoleQuad4QuarterSym.exo')
    max_p = [0., None]
    max_u = [0., None]
    for stage in F.stages.values():
        for increment in stage.increments:
            u = increment.field_outputs['U']
            for value in u.values:
                u1 = value.magnitude
                if max_u[0] < u1:
                    max_u = [u1, value]

            s = increment.field_outputs['S']
            for value in s.values:
                s1 = value.max_principal
                if max_p[0] < s1:
                    max_p = [s1, value]

    # External and internal element numbers
    xel = max_p[1].label
    x = F.get_elem_coord(xel)
    if plot:
        #print(max_p[0])
        V.Plot2D(deformed=1)

def test(plot=False):
    QuarterPlateDemo(plot=plot)

if __name__ == '__main__':
    test(plot=True)
