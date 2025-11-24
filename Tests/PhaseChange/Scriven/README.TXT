This folder includes main.cpp for Scriben bubble growth problem.
1. main.cpp
   Scriven problem: water at 1 bar, superheat = 1.25 K, initial radius = 50 microns
   If you comment out #define USE_VOF (i.e., //#define USE_VOF),
   then CIPCSL2 will be used. 
2. run
   $ mpiexec -np 8 ./Boil 120 < /dev/null > log.txt &
      input: dT1.25_R50microns.txt (temperature profile of initial condition)
      output: *.dat, *.bck
2. post-process
   $ post.sh
      input: log.txt, *.gnu
      output: radius.png
