#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from numpy import *
from pyfem2 import *

# READ MESH FROM FILE
mesh = GenesisMesh('PlateWithHoleTria3.g')

# ELASTIC MATERIAL MODEL
mat = Material('Material-1', elastic={'E':10e6, 'Nu':.29})

# CREATE THE MODEL AND ASSIGN PROPERTIES
V = FiniteElementModel(mesh=mesh, jobid='Plane2')
V.AssignProperties('ElementBlock1', PlaneStrainTria3, mat, t=1)

# PRESCRIBE FIXED BCS TO MODEL
V.PrescribedBC('LeftHandSide', X)
V.PrescribedBC('BottomLeft', Y)

# CREATE LOAD STEP AND PRESCRIBED NONHOMOGENEOUS BCS TO IT
step = V.StaticStep()
step.PrescribedBC('RightHandSide', X, .1)

# RUN THE STEP TO SOLVE FOR THE UNKNOWN DEGREES OF FREEDOM
step.run()

# WRITE RESULTS
V.WriteResults()

if not os.getenv('NOGRAPHICS'):
    # VISUALIZE RESULTS
    V.Plot2D(show=1, deformed=1)
