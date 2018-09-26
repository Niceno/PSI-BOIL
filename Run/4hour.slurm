#!/bin/bash
#SBATCH --ntasks=64
#SBATCH --time=03:59:00

# set some parameters
TEST=time.txt
OUT=out
END=120000
NP=64
TIME=240

# go to working directory
cd /scratch/rosa/bniceno/Boxed

# if TEST exists, read it and launch the job if needed
if [ -f $TEST ]; then
  read TIMESTEP < $TEST
  if [ "$TIMESTEP" -lt "$END" ]; then
    aprun -n $NP ./Boil >> $OUT
    sbatch $0
  fi

# if TEST file does not exist, launch the first job
else
  aprun -n $NP ./Boil > $OUT
  sbatch $0
fi

exit
