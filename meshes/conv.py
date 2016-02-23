import os
import glob
files = glob.glob('./meshes/*.vtu')
from pyfem2.mesh import VTU2Genesis

for filename in files:
    if os.path.isfile(os.path.splitext(filename)[0]+'.g'):
        continue
    print filename
    VTU2Genesis(filename=filename)
