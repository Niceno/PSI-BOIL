import sys
noise_file = "uvw-c_pall_*.plt database"
 
OpenDatabase(noise_file)

# Contour
AddPlot("Contour", "A")
c = ContourAttributes()
c.colorType = 0
c.contourMethod = 1
c.contourValue = 0.5
c.singleColor = (255, 255, 255, 255)
SetPlotOptions(c)

# View
v = GetView3D()
v.viewNormal = (0, -1, 0)
v.viewUp = (0, 0, 1)
v.perspective = 0
SetView3D(v)
DrawPlots()

for i in range(GetDatabaseNStates()):
  SaveWindow()
  TimeSliderNextState()

sys.exit() 
