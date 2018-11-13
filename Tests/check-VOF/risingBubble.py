import sys
noise_file = "uvw-c-press_pall_000004.plt"
 
OpenDatabase(noise_file)

# Contour
AddPlot("Contour", "A")
c = ContourAttributes()
c.colorType = 0
c.contourMethod = 1
c.contourValue = 0.5
c.singleColor = (255, 255, 255, 128)
SetPlotOptions(c)

# Vector
DefineVectorExpression("tension", "{U,V,W}")
AddPlot("Vector", "tension")
v = VectorAttributes()
v.useStride = 1
v.stride = 1
SetPlotOptions(v)
AddOperator("Slice")
s = SliceAttributes()
s.project2d = 0
SetOperatorOptions(s)

# Vector Module
DefineScalarExpression("Module","sqrt(U^2+V^2+W^2)")
AddPlot("Pseudocolor", "Module")
AddOperator("Slice")
s = SliceAttributes()
s.project2d = 0
SetOperatorOptions(s)

# View
v = GetView3D()
v.viewNormal = (0, -1, 0)
v.viewUp = (0, 0, 1)
v.perspective = 0
SetView3D(v)

w = SaveWindowAttributes()
w.format = w.PNG
w.width, w.height = 2048, 2048
SetSaveWindowAttributes(w)

DrawPlots()
SaveWindow()
sys.exit() 
