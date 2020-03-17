from numpy import sqrt
from felab.fe_model import FEModel
from felab.elemlib import DC2D3
from felab.constants import T, ALL, ILO, IHI, JLO, JHI
from felab.io.plot import plot2d


def demo_heat(plot=False):
    # Create the model
    V = FEModel(jobid="Heat1")

    # Read mesh from file
    V.genesis_mesh("./data/PlateWithHoleTria3Fine.g")

    # Create a material and define the thermal conductivity
    mat = V.material("Material-1")
    mat.isotropic_thermal_conductivity(12)

    # Define an alement block of diffusive heat transfer elements with material mat
    V.assign_properties(
        element_block="ElementBlock1", element_type=DC2D3, material=mat, t=1
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
