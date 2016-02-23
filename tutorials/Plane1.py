#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *
import matplotlib.pyplot as plt

V = Plane2DModel()
V.GenesisMesh('../meshes/QuarterCylinderQuad4.g')
V.Material('Material-1')
V.materials['Material-1'].Elastic(mu=1, Nu=.499)
V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)
V.PrescribedBC('Nodeset-200', X)
V.PrescribedBC('Nodeset-201', Y)
V.Pressure('Surface-1', 1.)
V.Solve()
V.WriteResults('VolumeLocking.Linear')
