#!/bin/bash

now=$(date +'%Y-%m-%d')
echo $now
rm -rf $now
mkdir $now
mkdir $now/Srcfiles
mkdir $now/Execs
mkdir $now/Srcfiles/VOF
mkdir $now/Srcfiles/PhaseChangeVOF
mkdir $now/Srcfiles/EnthalpyFD

cp ../Equation/Centered/VOF/*.cpp $now/Srcfiles/VOF/
cp ../Equation/Centered/VOF/*.h   $now/Srcfiles/VOF/
cp ../Equation/Centered/PhaseChangeVOF/*.cpp $now/Srcfiles/PhaseChangeVOF/
cp ../Equation/Centered/PhaseChangeVOF/*.h   $now/Srcfiles/PhaseChangeVOF/
cp ../Equation/Centered/EnthalpyFD/*.cpp   $now/Srcfiles/EnthalpyFD/
cp ../Equation/Centered/EnthalpyFD/*.h     $now/Srcfiles/EnthalpyFD/

cd ..

echo '##############'
echo '### Stefan ###'
echo '##############'
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

echo '###############'
echo '### Sucking ###'
echo '###############'
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

echo '####################'
echo '### BubbleGrowth ###'
echo '####################'
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
