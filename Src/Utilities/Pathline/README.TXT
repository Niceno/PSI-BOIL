Description:
pathline.exe import particle files (e.g. particles_000*.dat) and create pathlines for each particle.
input:  particle data written in Tecplot ASCII format (particles_000*.dat)
output: pathline.dat (Tecplot format)

1. compile
> ifort pathline.f90 -o pathline.exe

2. execute
> pathline.exe
