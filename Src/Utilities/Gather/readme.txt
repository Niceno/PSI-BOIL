tecio64.a  is created in tecio (old library using TecDat110)
libtecio.a is created in teciosrc (new library using TecDat142)

0. Compile libtecio.a
 cd teciosrc
 read "readme.txt" and compile.
 cmake .
 make

1. Compile gather.f90
  -DVISIT:  without this option, Visit cannot read plt files
  -DZIP:    gzip after the output of plt file
  -DTEC142: with this option, tecio version 142 is used. Without this option, tecio version 110 is used.
            tecio 142 is faster than 110. But plt-file output from 110 can be read directly 
  -DCSCS:   include memory.f90 which may be useful for debug
1.1 Intel compiler with OpenMP
 ifort -fpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -qopenmp -o gather.exe
 ifort -fpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -qopenmp -o gather-zip.exe
 ifort -fpp -DVISIT -DTEC142 gather.f90 ./libtecio.a -lm -lstdc++ -qopenmp -o gather-tec142.exe

1.2 GNU compiler with OpenMP
 gfortran -cpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -fopenmp -fcray-pointer -o gather.exe
 gfortran -cpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -fopenmp -fcray-pointer -o gather-zip.exe
 gfortran -cpp -DVISIT -DTEC142 gather.f90 ./libtecio.a -lm -lstdc++ -fopenmp -fcray-pointer -o gather-tec142.exe

1.3 On eiger.cscs.ch
 module switch PrgEnv-cray PrgEnv-intel
   ftn -fpp -DVISIT -DCSCS gather.f90 ./tecio64.a -lm -lstdc++ -qopenmp -o gather.exe memory.f90
   ftn -fpp -DVISIT -DTEC142 -DCSCS gather.f90 ./libtecio.a -lm -lstdc++ -qopenmp -o gather.exe memory.f90
 module switch PrgEnv-cray PrgEnv-gnu
   ftn -cpp -DVISIT -DCSCS gather.f90 ./tecio64.a -fcray-pointer -lm -lstdc++ -fopenmp -o gather.exe memory.f90
   ftn -cpp -DVISIT -DTEC142 -DCSCS gather.f90 ./libtecio.a -fcray-pointer -lm -lstdc++ -fopenmp -o gather-tec142.exe memory.f90

3. Executable files
   gather.exe: output binary file is tec110 format
   gather-zip.exe: output binary file is tec110 format + gzip
   gather-tec142.exe: output binary file is tec142 format 

4. Options for gather.exe
   --debug (= verbose, not slow down)
   --wait 10  (wait 10 min before stop)

5. Compile preplot
 g++ preplot.cpp -DPLOT3D -DUNIXX -DLINUX -o preplot
