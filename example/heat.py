import numpy as np

from felab import *
from felab.elemlib import DC2D3


def demo_plate_with_hole_fine(plot=False):
    mesh = genesis_mesh("./data/PlateWithHoleTria3Fine.g")
    V = FEModel(mesh=mesh)
    k, h, Too = 12, 250, 25
    fun = lambda x: 1000 / np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)
    mat = Material(name="Material-1", thermal_conductivity=k)
    V.assign_properties(element_block="ElementBlock1", element_type=DC2D3, material=mat)
    step = V.heat_transfer_step()
    step.dirichlet_bc(ILO, T, 200)
    step.dirichlet_bc(IHI, T, 50)
    step.dflux(JLO, 2000)
    step.sfilm(JHI, Too, h)
    step.HeatGeneration(ALL, fun)
    step.run()
    if plot:
        V.mesh.PlotScalar2D(step.dofs.flatten())


def demo_plate_with_hole_coarse():
    k, h, Too = 12, 250, 25  # noqa: F841
    mesh = genesis_mesh("./data/PlateWithHoleTria3.g")
    V = FEModel(jobid="HeatPlateWithHole", mesh=mesh)
    mat = Material(name="Material-1", thermal_conductivity=k)
    V.assign_properties(element_block="ElementBlock1", element_type=DC2D3, material=mat)
    V.initial_temperature(ALL, 50)
    step = V.heat_transfer_step()
    step.dirichlet_bc(BOUNDARY, T, 50)
    fun = lambda x: 1000 / np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)
    step.HeatGeneration(ALL, fun)
    step.run()
    V.write_results()


if __name__ == "__main__":
    demo_plate_with_hole_fine(plot=True)
    demo_plate_with_hole_coarse()
