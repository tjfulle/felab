import numpy as np

from felab.fe_model import fe_model
from felab.elemlib import DC2D3
from felab.constants import T, ALL, ILO, IHI, JLO, JHI, BOUNDARY


def demo_plate_with_hole_fine(plot=False):
    V = fe_model()
    V.genesis_mesh("./data/PlateWithHoleTria3Fine.g")
    k, h, Too = 12, 250, 25
    fun = lambda x: 1000 / np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(k)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(ILO, T, 200)
    step.assign_prescribed_bc(IHI, T, 50)
    step.SurfaceFlux(JLO, 2000)
    step.SurfaceConvection(JHI, Too, h)
    step.HeatGeneration(ALL, fun)
    step.run()
    if plot:
        V.mesh.PlotScalar2D(step.dofs.flatten())


def demo_plate_with_hole_coarse():
    k, h, Too = 12, 250, 25  # noqa: F841
    V = fe_model(jobid="HeatPlateWithHole")
    V.genesis_mesh("./data/PlateWithHoleTria3.g")
    V.create_material("Material-1")
    V.materials["Material-1"].isotropic_thermal_conductivity(k)
    V.assign_properties("ElementBlock1", DC2D3, "Material-1")
    V.assign_initial_temperature(ALL, 50)
    step = V.create_heat_transfer_step()
    step.assign_prescribed_bc(BOUNDARY, T, 50)
    fun = lambda x: 1000 / np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)
    step.HeatGeneration(ALL, fun)
    step.run()
    V.write_results()


if __name__ == "__main__":
    demo_plate_with_hole_fine(plot=True)
    demo_plate_with_hole_coarse()
