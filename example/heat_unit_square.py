from numpy import sqrt
from felab.fe_model import fe_model
from felab.elemlib import DC2D3
from felab.constants import ILO, IHI, JLO, T


def test_heat_unit_square(plot=False):
    # Create the model
    V = fe_model(jobid="Heat1")

    # Read mesh from file
    V.abaqus_mesh("./data/mesh.inp")

    # Create a material and define the thermal conductivity
    mat = V.create_material("Material-1")
    mat.isotropic_thermal_conductivity(12)
    print(V.mesh.element_blocks[0].elecon)

    # Define an alement block of diffusive heat transfer elements with material mat
    V.assign_properties("EALL", DC2D3, mat)

    # Fix temperatures on left and right edge
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(ILO, T, 200)
    step.assign_prescribed_bc(IHI, T, 50)

    # Define surface flux on bottome edge of domain
    step.SurfaceFlux(JLO, 2000)

    # Define surface convection on top edge of domain
    # Too, h = 25, 250
    # step.SurfaceConvection(JHI, Too, h)

    # Define a function specifying the heat generation
    def fun(x):
        return 1000.0 / sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)

    # step.HeatGeneration(ALL, fun)

    # Solve for the unknown degrees of freedom
    step.run()
    step.print_stiffness_structure()

    V.write_results()

    if plot:
        V.Plot2D(colorby=step.dofs.flatten(), show=1)

        # PlotScalar2D(V.mesh.coord, V.mesh.elecon, V.dofs.flatten())


if __name__ == "__main__":
    demo_heat_unit_square(plot=True)
