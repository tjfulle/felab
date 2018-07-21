from felab import *

def PostProcessDemo():
    mesh = genesis_mesh('./data/PlateWithHoleQuad4.g')
    mat = Material('Material-1')
    mat.elastic(E=100, Nu=.2)

    V = fe_model(mesh=mesh, jobid='PlateWithHoleQuad4')

    V.assign_properties('ElementBlock1', CPE4, mat, t=1)
    V.assign_prescribed_bc('LeftHandSide', X)
    V.fix_nodes('PinNode')

    step = V.create_static_step()
    step.assign_surface_load(IHI, [1,0])
    step.run()

    V.write_results()

    F = File('PlateWithHoleQuad4.exo')
    max_p = [0., None]
    max_u = [0., None]
    for step in F.steps.values():
        for frame in step.frames:
            u = frame.field_outputs['U']
            for value in u.values:
                u1 = value.magnitude
                if max_u[0] < u1:
                    max_u = [u1, value]

            s = frame.field_outputs['S']
            for value in s.values:
                s1 = value.max_principal
                if max_p[0] < s1:
                    max_p = [s1, value]

    # External and internal element numbers
    xel = max_p[1].label
    x = F.get_elem_coord(xel)


def test():
    PostProcessDemo()

if __name__ == '__main__':
    test()
