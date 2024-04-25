1. compile main.cpp
2. run (e.g. >mpiexec -np 16 ./Boil < /dev/null > log.txt &)
3. post-process
   > dogather
   > pathline.exe < input-pathline.txt
   > dopreplot
   > tecmcr.exe < input-tecmcr.txt
   NOTE: pathline.exe is in ./Src/Utilities/Pathline/pathline.f90
         tecmcr.exe is in ./Src/Utilities/TecplotMacro/tecmcr.f90
   > tec360 3D.lay &
     tec360 GUI: Scripting (top bar) > Play Macro/Script: select png.mcr

        
