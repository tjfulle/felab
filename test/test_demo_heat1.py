from numpy import *
from felab import *


def Heat1(plot=False):
    # Create the model
    V = fe_model(jobid="Heat1")

    # Read mesh from file
    V.genesis_mesh("./data/PlateWithHoleTria3Fine.g")

    # Create a material and define the thermal conductivity
    mat = V.create_material("Material-1")
    mat.isotropic_thermal_conductivity(12)

    # Define an alement block of diffusive heat transfer elements with material mat
    V.assign_properties("ElementBlock1", DC2D3, mat)

    # Fix temperatures on left and right edge
    stage = V.create_heat_transfer_stage()
    stage.assign_prescribed_bc(ILO, T, 200)
    stage.assign_prescribed_bc(IHI, T, 50)

    # Define surface flux on bottome edge of domain
    stage.SurfaceFlux(JLO, 2000)

    # Define surface convection on top edge of domain
    Too, h = 25, 250
    stage.SurfaceConvection(JHI, Too, h)

    # Define a function specifying the heat generation
    def fun(x):
        return 1000.0 / sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)

    stage.HeatGeneration(ALL, fun)

    # Solve for the unknown degrees of freedom
    stage.run()

    V.write_results()

    if plot:
        V.Plot2D(colorby=stage.dofs.flatten(), show=1)
        # PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())


def test(plot=False):
    Heat1(plot=plot)


if __name__ == "__main__":
    test(plot=True)
