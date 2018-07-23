from felab import *

mu = 1.
Nu = .499
E = 2. * mu * (1. + Nu)

def LinearSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.Linear')
    V.abaqus_mesh('./data/ThickCylinder_Linear.inp')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ALL', CPE4R, mat, t=1)
    V.assign_prescribed_bc('SymYZ', X)
    V.assign_prescribed_bc('SymXZ', Y)
    stage = V.create_static_stage()
    # Pressure on inside face
    stage.assign_pressure('SurfID', 1.)
    stage.run()
    V.write_results()
    if ax is not None:
        ax = V.Plot2D(deformed=1, color='b', linestyle='-.', ax=ax, label='Linear',
                      xlim=(-.2,5), ylim=(-.2, 5))
        return ax
    return None

def QuadraticSolution(ax=None):
    V = fe_model(jobid='VolumeLocking.Quadratic')
    V.abaqus_mesh('./data/ThickCylinder_Quadratic.inp')
    mat = V.create_material('Material-1')
    mat.elastic(E=E, Nu=Nu)
    V.assign_properties('ALL', CPE8B, mat, t=1)
    V.assign_prescribed_bc('SymYZ', X)
    V.assign_prescribed_bc('SymXZ', Y)
    # Pressure on inside face
    stage = V.create_static_stage()
    stage.assign_pressure('SurfID', 1.)
    stage.run()
    V.write_results()
    if ax is not None:
        return ax
    return None

def WriteAnalyticSolution(plot=False):
    mesh = Mesh(filename='./data/ThickCylinder_Linear.inp')
    ix = where(mesh.coord[:,1]<=1e-12)
    a = mesh.coord[ix][:,0].min()
    b = mesh.coord[ix][:,0].max()
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
        ax = mesh.Plot2D(xy=u+mesh.coord, ax=None, color='orange', weight=6,
                         label='Analytic')
        return ax
    return None

def test(plot=False):
    ax = None
    ax = WriteAnalyticSolution(plot=plot)
    ax = LinearSolution(ax)
    QuadraticSolution()

    if plot:
        import matplotlib.pyplot as plt
        plt.legend()
        plt.show()


if __name__ == '__main__':
    test(plot=True)
