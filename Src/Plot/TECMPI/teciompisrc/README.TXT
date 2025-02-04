on merlin and CSCS
$ module load Boost
$ cmake .
$ make

on euler
$ module load boost
$ cmake .
$ make

on merlin7
$ module load cmake
$ module list
Currently Loaded Modules:
 1) gcc14.2.0           2) cuda12.2.0          3) libfabric1.15.2.0  
 4) libevent2.1.12      5) pmix5.0.3           6) openmpi5.0.5 
$ cmake .
$ make

re-compile
$ make clean
$ rm CMakeCache.txt
$ rm -rf CMakeFiles/ 
