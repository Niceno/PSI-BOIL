import sys
noise_file = "uvw-c_pall_*.plt database"
 
OpenDatabase(noise_file)

AddPlot("Pseudocolor", "A")
DrawPlots()

for i in range(GetDatabaseNStates()):
  SaveWindow()
  TimeSliderNextState()
 
sys.exit()
