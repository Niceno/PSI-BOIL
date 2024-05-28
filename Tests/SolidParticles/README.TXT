1. compile main.cpp
2. run (e.g. >mpiexec -np 16 ./Boil < /dev/null > log.txt &)
3. post-process: file conversion
   > dogather
   > pathline.exe < input-pathline.txt
   > dopreplot
4A. post-process: visualization using Tecplot
   > tecmcr.exe < input-tecmcr.txt
   NOTE: pathline.exe is in ./Src/Utilities/Pathline/pathline.f90
         tecmcr.exe is in ./Src/Utilities/TecplotMacro/tecmcr.f90
   > tec360 3D.lay &
     tec360 GUI: Scripting (top bar) > Play Macro/Script: select png.mcr
4B. post-process: visualization usinv visit
   > visit
     Open: uvw*.plt, particles_*.dat, pathlines.dat
     (Warning! When you open particles_*dat files, specify the file type as Tecplot)
     particles can be visualized by Add > mesh > points > mesh
     pathlines can be visualized by Add > mesh
        
