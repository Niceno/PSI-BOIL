0. Compile libtecio.a
 cd teciosrc
 read "readme.txt" and compile.
 cmake .
 make

1. Compile gather.f90
 ifort -fpp -DVISIT gather.f90 ./libtecio.a -lm -lstdc++ -o gather.exe
 ifort -fpp -DVISIT -DZIP gather.f90 ./libtecio.a -lm -lstdc++ -o gather-zip.exe
 ifort -fpp -DVISIT -DSZPLT gather.f90 ./libtecio.a -lm -lstdc++ -o gather-szplt.exe

2. Compile preplot
 g++ preplot.cpp -DPLOT3D -DUNIXX -DLINUX -o preplot

3. Compile gather.f90 on ROSA
 module switch PrgEnv-gnu PrgEnv-intel
 ftn gather.f90 ./libtecio.a -fpp -lm -lstdc++ -o gather.exe
 ftn gather.f90 ./libtecio.a -fpp -DZIP -lm -lstdc++ -o gather-zip.exe
 
 use module PrgEnv-pgi
 ftn gather.f90 ./libtecio.a -cpp -lm -lstdc++ -fcray-pointer -o gather.exe
 ftn gather.f90 ./libtecio.a -cpp -DZIP -lm -lstdc++ -fcray-pointer -o gather-zip.exe
