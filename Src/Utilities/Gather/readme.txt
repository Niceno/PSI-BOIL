tecio64.a  is created in tecio (old library using TecDat110)
libtecio.a is created in teciosrc (new library using TecDat142)

0. Compile libtecio.a
 cd teciosrc
 read "readme.txt" and compile.
 cmake .
 make

1. Compile gather.f90
1.1 Intel compiler with OpenMP
 ifort -fpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -qopenmp -o gather.exe
 ifort -fpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -qopenmp -o gather-zip.exe
 ifort -fpp -DVISIT -DSZPLT gather.f90 ./libtecio.a -lm -lstdc++ -qopenmp -o gather-szplt.exe
1.2 GNU compiler with OpenMP
 gfortran -cpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -fopenmp -fcray-pointer -o gather.exe
 gfortran -cpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -fopenmp -fcray-pointer -o gather-zip.exe
 gfortran -cpp -DVISIT -DSZPLT gather.f90 ./libtecio.a -lm -lstdc++ -fopenmp -fcray-pointer -o gather-szplt.exe
1.3 On eiger.cscs.ch
 module switch PrgEnv-cray PrgEnv-intel
   ftn -fpp -DVISIT -DCSCS gather.f90 ./tecio64.a -lm -lstdc++ -qopenmp -o gather.exe memory.f90
 module switch PrgEnv-cray PrgEnv-gnu
   ftn -cpp -DVISIT -DCSCS gather.f90 ./tecio64.a -fcray-pointer -lm -lstdc++ -fopenmp -o gather.exe memory.f90

2. Execute gather.exe
 possible option --debug (= verbose, not slow down)

3. Compile preplot
 g++ preplot.cpp -DPLOT3D -DUNIXX -DLINUX -o preplot
