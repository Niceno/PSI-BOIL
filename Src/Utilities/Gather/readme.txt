tecio64.a  is created in tecio (old library using TecDat110)
libtecio.a is created in teciosrc (new library using TecDat142)

0. Compile libtecio.a
 cd teciosrc
 read "readme.txt" and compile.
 cmake .
 make

1. Compile gather.f90
1.1 Intel compiler
 ifort -fpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -o gather.exe
 ifort -fpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -o gather-zip.exe
 ifort -fpp -DVISIT -DSZPLT gather.f90 ./libtecio.a -lm -lstdc++ -o gather-szplt.exe
1.2 GNU compiler
 gfortran -cpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -fcray-pointer -o gather.exe
 gfortran -cpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -fcray-pointer -o gather-zip.exe
 gfortran -cpp -DVISIT -DSZPLT gather.f90 ./libtecio.a -lm -lstdc++ -fcray-pointer -o gather-szplt.exe

2. Compile preplot
 g++ preplot.cpp -DPLOT3D -DUNIXX -DLINUX -o preplot

3. Compile gather.f90 on ROSA
 module switch PrgEnv-gnu PrgEnv-intel
 ftn gather.f90 ./tecio64.a -fpp -lm -lstdc++ -o gather.exe
 
 use module PrgEnv-pgi
 ftn gather.f90 ./tecio64.a -cpp -lm -lstdc++ -fcray-pointer -o gather.exe
