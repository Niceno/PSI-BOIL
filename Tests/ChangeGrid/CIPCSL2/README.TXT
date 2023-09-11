1. compile original grid (main-org.cpp)
   run 20 time steps
   output: *20.bck, g1.grid
2. compile new grid (main-new.cpp)
   run up to time loop
   output: g2.grid
3. ChangeGrid
   input: g1.grid, g2.grid, *20.bck, input.txt
   output: *20.bck2
4. change file extension from *bck2 to *bck
   Run Boil for new grid

