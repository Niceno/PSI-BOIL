import sys
noise_file = "uvw-c-tpr-mdot_002000.plt"
 
OpenDatabase(noise_file)

# Pseudocolor
AddPlot("Pseudocolor", "A")

# View
v = GetView3D()
v.viewNormal = (1, 1, 1)
v.viewUp = (0, 0, 1)
v.perspective = 1
SetView3D(v)

DrawPlots()
SaveWindow()
sys.exit() 
