import sys
noise_file = "uvw-c-press_pall_000858.plt"
 
OpenDatabase(noise_file)

#Pseudocolor
AddPlot("Pseudocolor", "A")
AddOperator("Slice")
a = SliceAttributes()
a.originType = a.Intercept
a.originIntercept = 0.00098
SetOperatorOptions(a)

#Vector
DefineVectorExpression("tension", "{U,V,W}")
AddPlot("Vector", "tension")
AddOperator("Slice")
a = SliceAttributes()
a.originType = a.Intercept
a.originIntercept = 0.00098
SetOperatorOptions(a)

v = VectorAttributes()
v.useStride = 1
v.stride = 1
SetPlotOptions(v)

DrawPlots()
SaveWindow()
sys.exit() 
