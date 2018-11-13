import sys
noise_file = "uvw-c_pall_000400.plt"
 
OpenDatabase(noise_file)

AddPlot("Pseudocolor", "A")
AddOperator("Slice")
a = SliceAttributes()
a.originType = a.Intercept
a.originIntercept = 0.1
SetOperatorOptions(a)

DrawPlots()
SaveWindow()
sys.exit() 
