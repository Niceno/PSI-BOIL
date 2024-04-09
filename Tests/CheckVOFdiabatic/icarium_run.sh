#!/bin/bash

if [[ $1 == "" ]]; then
  echo "Path to running folder is required!"
  exit
fi

rootf=$PWD
echo $rootf
now=$(date +'%Y-%m-%d')
echo $now
rm -rf $now
mkdir $now
mkdir $now/Srcfiles
mkdir $now/Execs
mkdir $now/Srcfiles/VOF
mkdir $now/Srcfiles/PhaseChangeVOF
mkdir $now/Srcfiles/EnthalpyTif

cp ../Equation/Centered/VOF/*.cpp $now/Srcfiles/VOF/
cp ../Equation/Centered/VOF/*.h   $now/Srcfiles/VOF/
cp ../Equation/Centered/PhaseChangeVOF/*.cpp $now/Srcfiles/PhaseChangeVOF/
cp ../Equation/Centered/PhaseChangeVOF/*.h   $now/Srcfiles/PhaseChangeVOF/
cp ../Equation/Centered/EnthalpyTif/*.cpp   $now/Srcfiles/EnthalpyTif/
cp ../Equation/Centered/EnthalpyTif/*.h     $now/Srcfiles/EnthalpyTif/

cd ..

echo '############## Compilation starts ##############'

echo '### Stefan ###'
mkdir CheckVOFdiabatic/$now/Execs/Stefan
cp CheckVOFdiabatic/Mainfiles/main.cpp.stefan.lvl2 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/Stefan/Boil.lvl2
cp CheckVOFdiabatic/Mainfiles/main.cpp.stefan.lvl4 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/Stefan/Boil.lvl4
cp CheckVOFdiabatic/Mainfiles/main.cpp.stefan.lvl6 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/Stefan/Boil.lvl6

echo '### Sucking ###'
mkdir CheckVOFdiabatic/$now/Execs/Sucking
cp CheckVOFdiabatic/Mainfiles/main.cpp.sucking.lvl2 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/Sucking/Boil.lvl2
cp CheckVOFdiabatic/Mainfiles/main.cpp.sucking.lvl4 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/Sucking/Boil.lvl4
cp CheckVOFdiabatic/Mainfiles/main.cpp.sucking.lvl6 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/Sucking/Boil.lvl6

echo '### BubbleGrowth ###'
cp Equation/Staggered/Momentum/momentum_scale_outlet_velocity.cpp Equation/Staggered/Momentum/momentum_scale_outlet_velocity.cpp.old
cp CheckVOFdiabatic/Mainfiles/momentum_scale_outlet_velocity.cpp.spherical Equation/Staggered/Momentum/momentum_scale_outlet_velocity.cpp

mkdir CheckVOFdiabatic/$now/Execs/BubbleGrowth
cp CheckVOFdiabatic/Mainfiles/main.cpp.bubblegrowth.lvl2 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/BubbleGrowth/Boil.lvl2
cp CheckVOFdiabatic/Mainfiles/main.cpp.bubblegrowth.lvl4 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/BubbleGrowth/Boil.lvl4
cp CheckVOFdiabatic/Mainfiles/main.cpp.bubblegrowth.lvl6 main.cpp
make >& /dev/null
cp Boil CheckVOFdiabatic/$now/Execs/BubbleGrowth/Boil.lvl6

mv Equation/Staggered/Momentum/momentum_scale_outlet_velocity.cpp.old Equation/Staggered/Momentum/momentum_scale_outlet_velocity.cpp

echo '############## Compilation finished ##############'
cd CheckVOFdiabatic/$now/Execs
execloc=$PWD
echo $1
cd $1
echo $execloc
rm -rf $now
mkdir $now
cd $now

echo '############## Copying starts ##############'

echo '### Stefan ###'
mkdir 'Stefan'
cd 'Stefan'
cp $rootf/Scripts/stefan.py .
echo '#lvl2#'
mkdir 'lvl2'
cd 'lvl2'
cp $execloc/Stefan/Boil.lvl2 Boil
cd ..
echo '#lvl4#'
mkdir 'lvl4'
cd 'lvl4'
cp $execloc/Stefan/Boil.lvl4 Boil
cd ..
echo '#lvl6#'
mkdir 'lvl6'
cd 'lvl6'
cp $execloc/Stefan/Boil.lvl6 Boil
cd ../../

echo '### Sucking ###'
mkdir 'Sucking'
cd 'Sucking'
cp $rootf/Scripts/sucking.py .
cp $rootf/Scripts/sucking.theoretical.1 theor.new
cp $rootf/Scripts/sucking.theoretical.2 theor.err
echo '#lvl2#'
mkdir 'lvl2'
cd 'lvl2'
cp $execloc/Sucking/Boil.lvl2 Boil
cp $rootf/Mainfiles/input.sucking input.txt
cd ..
echo '#lvl4#'
mkdir 'lvl4'
cd 'lvl4'
cp $execloc/Sucking/Boil.lvl4 Boil
cp $rootf/Mainfiles/input.sucking input.txt
cd ..
echo '#lvl6#'
mkdir 'lvl6'
cd 'lvl6'
cp $execloc/Sucking/Boil.lvl6 Boil
cp $rootf/Mainfiles/input.sucking input.txt
cd ../../

echo '### BubbleGrowth ###'
mkdir 'BubbleGrowth'
cd 'BubbleGrowth'
cp $rootf/Scripts/bubblegrowth.betag.py betag.py
cp $rootf/Scripts/bubblegrowth.py .
echo '#lvl2#'
mkdir 'lvl2'
cd 'lvl2'
cp $execloc/BubbleGrowth/Boil.lvl2 Boil
cd ..
echo '#lvl4#'
mkdir 'lvl4'
cd 'lvl4'
cp $execloc/BubbleGrowth/Boil.lvl4 Boil
cd ..
echo '#lvl6#'
mkdir 'lvl6'
cd 'lvl6'
cp $execloc/BubbleGrowth/Boil.lvl6 Boil

echo '############## Copying ends ##############'
cd '../../'
echo $PWD
echo '############## Simulations start ##############'

echo '### Stefan ###'
cd 'Stefan'
echo '#lvl2#'
cd 'lvl2'
mpiexec --use-hwthread-cpus -np 1 ./Boil 10000 > log.txt
grep x-m log.txt > xmn.txt
for i in *dat; do preplot $i; done; rm *dat
echo '#lvl4#'
cd '../lvl4'
mpiexec --use-hwthread-cpus -np 1 ./Boil 10000 > log.txt
grep x-m log.txt > xmn.txt
for i in *dat; do preplot $i; done; rm *dat
echo '#lvl6#'
cd '../lvl6'
mpiexec --use-hwthread-cpus -np 1 ./Boil 10000 > log.txt
grep x-m log.txt > xmn.txt
for i in *dat; do preplot $i; done; rm *dat
cd ..
python3 stefan.py

echo '### Sucking ###'
cd '../Sucking'
echo '#lvl2#'
cd 'lvl2'
mpiexec --use-hwthread-cpus -np 1 ./Boil 10000 > log.txt
grep x-m log.txt > xmn.txt
for i in *dat; do preplot $i; done; rm *dat
echo '#lvl4#'
cd '../lvl4'
mpiexec --use-hwthread-cpus -np 1 ./Boil 10000 > log.txt
grep x-m log.txt > xmn.txt
for i in *dat; do preplot $i; done; rm *dat
echo '#lvl6#'
cd '../lvl6'
mpiexec --use-hwthread-cpus -np 1 ./Boil 10000 > log.txt
grep x-m log.txt > xmn.txt
for i in *dat; do preplot $i; done; rm *dat
cd ..
python3 sucking.py

echo '### BubbleGrowth ###'
cd '../BubbleGrowth'
echo '#lvl2#'
cd 'lvl2'
mpiexec --use-hwthread-cpus -np 32 ./Boil 10000 > log.txt
grep totalvol log.txt > tot.txt
dogather
echo '#lvl4#'
cd '../lvl4'
mpiexec --use-hwthread-cpus -np 32 ./Boil 10000 > log.txt
grep totalvol log.txt > tot.txt
dogather
echo '#lvl6#'
cd '../lvl6'
mpiexec --use-hwthread-cpus -np 32 ./Boil 10000 > log.txt
grep totalvol log.txt > tot.txt
dogather
cd ..
python3 bubblegrowth.py

echo '############## Simulations end ##############'


