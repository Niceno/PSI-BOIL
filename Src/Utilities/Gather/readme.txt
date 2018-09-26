1. Compile gather.f90
 ifort -fpp gather.f90 ./tecio64.a -lm -lstdc++ -o gather.exe
 ifort -fpp -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -o gather-zip.exe
 ifort -fpp -DVISIT gather.f90 ./tecio64.a -lm -lstdc++ -o gather-visit.exe
 ifort -fpp -DVISIT -DZIP gather.f90 ./tecio64.a -lm -lstdc++ -o gather-visit-zip.exe

2. Compile preplot
 g++ preplot.cpp -DPLOT3D -DUNIXX -DLINUX -o preplot

3. Compile gather.f90 on ROSA
 module switch PrgEnv-gnu PrgEnv-intel
 ftn gather.f90 ./tecio64.a -fpp -lm -lstdc++ -o gather.exe
 ftn gather.f90 ./tecio64.a -fpp -DZIP -lm -lstdc++ -o gather-zip.exe
 
 use module PrgEnv-pgi
 ftn gather.f90 ./tecio64.a -cpp -lm -lstdc++ -fcray-pointer -o gather.exe
 ftn gather.f90 ./tecio64.a -cpp -DZIP -lm -lstdc++ -fcray-pointer -o gather-zip.exe
