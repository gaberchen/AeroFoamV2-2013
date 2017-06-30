#!/usr/bin/python

# Import C++/Python wrapper
from AeroFoam import *

# Set path and case names
path = "/home/giulio/OpenFOAM/giulio-1.7.1/benchmarks/Airfoils"
case = "EdgeComparison"

# Initialization
Time    = pyOpenFOAM( path, case )
Mesh    = pyMesh()
Aero    = pySolver()
Elastic = pyInterface()
Plugin  = pyPlugin()

# Time loop
while Time.next():
    Mesh.next()
    Aero.next()
    Elastic.next()
    Plugin.next()
    
# End
pyEnd()    
