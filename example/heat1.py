from numpy import sqrt

from felab import *
from felab.io.plot import plot2d


def demo_heat(plot=False):
    # Read mesh from file
    mesh = genesis_mesh("./data/PlateWithHoleTria3Fine.g")

    # Create the model
    V = FEModel(jobid="Heat1", mesh=mesh)

    # Create a material and define the thermal conductivity
    mat = Material(name="Material-1", thermal_conductivity=12)

    # Define an alement block of diffusive heat transfer elements with material mat
    el = Element(type="DC2D3")
    V.assign_properties(
        element_block="ElementBlock1", element_type=el, material=mat, t=1
    )

    # Fix temperatures on left and right edge
    step = V.heat_transfer_step()
    step.dirichlet_bc(ILO, T, 200)
    step.dirichlet_bc(IHI, T, 50)

    # Define surface flux on bottome edge of domain
    step.dflux(JLO, 2000)

    # Define surface convection on top edge of domain
    Too, h = 25, 250
    step.sfilm(JHI, Too, h)

    # Define a function specifying the heat generation
    def fun(x):
        return 1000.0 / sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)

    step.HeatGeneration(ALL, fun)

    # Solve for the unknown degrees of freedom
    step.run()

    V.write_results()

    if plot:
        plot2d(model=V, colorby=step.dofs.flatten(), show=1)
        # PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())


if __name__ == "__main__":
    demo_heat(plot=True)
