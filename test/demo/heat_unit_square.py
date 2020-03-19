import numpy as np
from felab import *
from felab.io.plot import plot2d


def demo_heat_unit_square(data_path, plot=False):
    # Read mesh from file
    mesh = abaqus_mesh(os.path.join(data_path, "mesh.inp"))

    # Create the model
    V = FEModel(jobid="Heat1", mesh=mesh)

    # Create a material and define the thermal conductivity
    mat = Material(name="Material-1", thermal_conductivity=12)
    print(V.mesh.element_blocks[0].elecon)

    # Define an alement block of diffusive heat transfer elements with material mat
    el = Element(type="DC2D3")
    V.assign_properties(element_block="EALL", element_type=el, material=mat)

    # Fix temperatures on left and right edge
    step = V.heat_transfer_step()
    step.dirichlet_bc(ILO, T, 200)
    step.dirichlet_bc(IHI, T, 50)

    # Define surface flux on bottome edge of domain
    step.dflux(JLO, 2000)

    # Define surface convection on top edge of domain
    # Too, h = 25, 250
    # step.sfilm(JHI, Too, h)

    # Define a function specifying the heat generation
    def fun(x):
        return 1000.0 / np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)

    # step.HeatGeneration(ALL, fun)

    # Solve for the unknown degrees of freedom
    step.run()
    step.print_stiffness_structure()

    V.write_results()

    if plot:
        plot2d(model=V, colorby=step.dofs.flatten(), show=1)

        # PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())


if __name__ == "__main__":
    import os
    this_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(this_path, "..", "data")
    demo_heat_unit_square(data_path, plot=True)
