from felab import *

mu = 1.
Nu = .499
E = 2. * mu * (1. + Nu)

def LinearSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.Linear')
    V.genesis_mesh('./data/QuarterCylinderQuad4.g')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ElementBlock1', CPE4, mat, t=1)
    V.assign_prescribed_bc('Nodeset-200', X)
    V.assign_prescribed_bc('Nodeset-201', Y)
    # Pressure on inside face
    step = V.create_static_step()
    step.assign_pressure('SURFACE-1', 1.)
    step.run()
    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, color='b', linestyle='-.', ax=ax, label='Linear',
                      xlim=(-.2,5), ylim=(-.2, 5))
        return ax
    return None

def BBarSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.BBar')
    V.genesis_mesh('./data/QuarterCylinderQuad4.g')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ElementBlock1', CPE4B, mat, t=1)
    V.assign_prescribed_bc('Nodeset-200', X)
    V.assign_prescribed_bc('Nodeset-201', Y)
    # Pressure on inside face
    step = V.create_static_step()
    step.assign_pressure('SURFACE-1', 1.)
    step.run()
    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, ax=ax, color='b', label='bbar')
        return ax
    return None

def ReducedIntegrationSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.Reduced')
    V.genesis_mesh('./data/QuarterCylinderQuad4.g')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ElementBlock1', CPE4R, mat, t=1)
    V.assign_prescribed_bc('Nodeset-200', X)
    V.assign_prescribed_bc('Nodeset-201', Y)
    # Pressure on inside face
    step = V.create_static_step()
    step.assign_pressure('SURFACE-1', 1.)
    step.run()
    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, ax=ax, color='b', label='reduced')
        return ax
    return None

def SelReducedIntegrationSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.SelReduced')
    V.genesis_mesh('./data/QuarterCylinderQuad4.g')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ElementBlock1', CPE4RS,
                       mat, t=1)
    V.assign_prescribed_bc('Nodeset-200', X)
    V.assign_prescribed_bc('Nodeset-201', Y)
    # Pressure on inside face
    step = V.create_static_step()
    step.assign_pressure('SURFACE-1', 1.)
    step.run()
    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, ax=ax, color='g', label='Sel. reduced integration')
        return ax
    return None

def QuadraticSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.Quadratic')
    V.genesis_mesh('./data/QuarterCylinderQuad8.g')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ElementBlock1', CPE8B, mat, t=1)
    V.assign_prescribed_bc('Nodeset-200', X)
    V.assign_prescribed_bc('Nodeset-201', Y)
    # Pressure on inside face
    step = V.create_static_step(solver=NEWTON)
    step.assign_pressure('SURFACE-1', 1.)
    #V.assign_surface_load("SURFACE-300", [0.195090322, 0.98078528])
    #V.assign_surface_load("SURFACE-301", [0.555570233, 0.831469612])
    #V.assign_surface_load("SURFACE-302", [0.831469612, 0.555570233])
    #V.assign_surface_load("SURFACE-303", [0.98078528, 0.195090322])
    step.run()
    V.write_results()

def WriteAnalyticSolution(plot=False):
    mesh = Mesh(filename='./data/QuarterCylinderQuad4.g')
    a = mesh.coord[0, 1]
    b = mesh.coord[-1, 0]
    p = 1.
    u = zeros_like(mesh.coord)
    for (i, x) in enumerate(mesh.coord):
        r = sqrt(x[0] ** 2 + x[1] ** 2)
        term1 = (1. + Nu) * a ** 2 * b ** 2 * p / (E * (b ** 2 - a ** 2))
        term2 = 1. / r + (1. - 2. * Nu) * r / b ** 2
        ur = term1 * term2
        u[i, :] = ur * x[:] / r
    mesh.put_nodal_solution('VolumeLocking.Analytic', u)
    if plot:
        ax = mesh.Plot2D(xy=u+mesh.coord, ax=None, color='orange', weight=8,
                         label='Analytic')
        return ax
    return None


def test(plot=False):
    ax = None
    ax = WriteAnalyticSolution(plot=plot)
    ax = ReducedIntegrationSolution(ax)
    ax = SelReducedIntegrationSolution(ax)
    ax = LinearSolution(ax)
    QuadraticSolution()

    if plot:
        import matplotlib.pyplot as plt
        plt.legend()
        plt.show()


if __name__ == '__main__':
    test(plot=True)
