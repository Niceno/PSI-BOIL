1. Copy main.cpp and update_step.cpp to Src, and compile Boil.
2. Copy floor.stl and Boil to the directory for the calculation.
3. Execute Boil. 
   > mpiexec -np 8 ./Boil 1440 < /dev/null > log.txt
4. Post-process: Copy post.sh and *.gnu to the calculation directory
   > post.sh
   Three figures are generated: diameter.png, hflux.png, t_nucl.png

