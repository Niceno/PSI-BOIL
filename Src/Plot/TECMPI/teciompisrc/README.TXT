$ module load Boost
$ cmake .
$ make

on euler
$ module load boost
$ cmake .
$ make

recompile
$ make clean
$ rm CMakeCache.txt
$ rm -rf CMakeFiles/
$ make
