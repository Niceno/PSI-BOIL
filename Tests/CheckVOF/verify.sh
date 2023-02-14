#!/bin/tcsh

echo '### How to execute ################'
echo '%./verify.sh tec     > use tecplot for visualization'
echo '%./verify.sh visit   > use visit for visualization'
echo '%./verify.sh visit sbatch  > use sbatch for heavy jobs, visit for visualization'
echo '### End: how to execute ###########'

if ( $1 =~ "" || ($1 !~ "tec" && $1 !~ "visit") ) then
  echo "Wrong syntax!"
  echo "./verify.sh tec            > use tecplot for visualization"
  echo "./verify.sh visit          > use visit for visualization"
  echo '%./verify.sh visit sbatch  > use sbatch for heavy jobs, visit for visualization'
  exit
endif

echo $1
echo $2


### Files requred: *.txt *.gth *.stl, main*.cpp *.gnu *.lay

rm -rf Result_tmp
mkdir Result_tmp
mkdir Result_tmp/SrcVOF
cp ../Equation/Centered/VOF/*.cpp Result_tmp/SrcVOF
cp ../Equation/Centered/VOF/*.h   Result_tmp/SrcVOF
cd ..

echo '###############'
echo '### Zalesak ###'
echo '###############'
mkdir CheckVOF/Result_tmp/Zalesak
cp CheckVOF/main-zalesak.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Zalesak
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpiexec -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null

cp ../../zalesak.gth .
dogather
  
if ( $1 =~ "tec" ) then
  # tecplot     
  cp ../../make_png.mcr .
  cp ../../zalesak.lay .
  foreach i (*.plt)
    tec360 zalesak.lay -b make_png.mcr $i >& /dev/null
    mv tmp.png $i:r.png
  end
else
  # Visit
  cp ../../zalesak.py .
  visit -cli -nowin -s zalesak.py >& /dev/null
endif
cd ../../..

echo '##W##################'
echo '### Circle vortex ###'
echo '#####################'
mkdir CheckVOF/Result_tmp/CircleVortex
cp CheckVOF/main-circle-vortex.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/CircleVortex
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpiexec -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null

cp ../../circleVortex.gth .
dogather

if ( $1 =~ "tec" ) then
  # tecplot
  cp ../../make_png.mcr .
  cp ../../circleVortex.lay .
  foreach i (*.plt)
    tec360 circleVortex.lay -b make_png.mcr $i >& /dev/null
    mv tmp.png $i:r.png
  end
else
  # Visit
  cp ../../circleVortex.py .
  visit -cli -nowin -s circleVortex.py >& /dev/null
endif
cd ../../..

echo '###################'
echo '### Corner flow ###'
echo '###################'
mkdir CheckVOF/Result_tmp/CornerFlow
cp CheckVOF/main-corner-flow.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/CornerFlow
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpiexec -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null

cp ../../cornerFlow.gth .
dogather

if ( $1 =~ "tec" ) then
  # tecplot
  cp ../../make_png.mcr .
  cp ../../cornerFlow.lay .
  foreach i (*.plt)
    tec360 cornerFlow.lay -b make_png.mcr $i >& /dev/null
    mv tmp.png $i:r.png
  end
else
  # Visit
  cp ../../cornerFlow.py .
  visit -cli -nowin -s cornerFlow.py >& /dev/null
endif
cd ../../..

echo '################################'
echo '### sliding sphare near wall ###'
echo '################################'
mkdir CheckVOF/Result_tmp/Sliding_Sphare_nearWall
cp CheckVOF/main-slide-sphare-nearWall.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Sliding_Sphare_nearWall
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpiexec -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null

cp ../../slidingSphare.gth .
dogather

if ( $1 =~ "tec" ) then
  # tecplot 1
  cp ../../make_png.mcr .
  cp ../../slidingSphareWall-cont.lay .
  tec360 slidingSphareWall-cont.lay -b make_png.mcr >& /dev/null
  mv tmp.png slidingSphareWall-cont.png
else  
  # Visit 1
  cp ../../slidingSphareWall-cont.py .
  visit -cli -nowin -s slidingSphareWall-cont.py >& /dev/null
  mv visit0000.png slidingSphareWall-cont.png 
endif

if ( $1 =~ "tec" ) then 
  # tecplot 2
  cp ../../slidingSphareWall-iso.lay .
  tec360 slidingSphareWall-iso.lay -b make_png.mcr >& /dev/null
  mv tmp.png slidingSphareWall-iso.png
else
  # Visit 2
  cp ../../slidingSphareWall-iso.py .
  visit -cli -nowin -s slidingSphareWall-iso.py >& /dev/null
  mv visit0000.png slidingSphareWall-iso.png
endif
cd ../../..

echo '##############################'
echo '### sliding sphare near IB ###'
echo '##############################'
mkdir CheckVOF/Result_tmp/Sliding_Sphare_nearIB
cp CheckVOF/main-slide-sphare-nearIB.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Sliding_Sphare_nearIB
cp ../../../Boil .
cp ../../../main.cpp .
cp ../../floor.stl .
set sec0 = `date +%s`
mpiexec -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
cp ../../slidingSphare.gth .
dogather

if ( $1 =~ "tec" ) then
  # tecplot
  cp ../../make_png.mcr .
  cp ../../slidingSphareIB-cont.lay .
  tec360 slidingSphareIB-cont.lay -b make_png.mcr >& /dev/null
  mv tmp.png slidingSphareIB-cont.png
  cp ../../slidingSphareIB-iso.lay .
  tec360 slidingSphareIB-iso.lay -b make_png.mcr >& /dev/null
  mv tmp.png slidingSphareIB-iso.png
else
  # Visit
  cp ../../slidingSphareIB-cont.py .
  visit -cli -nowin -s slidingSphareIB-cont.py >& /dev/null
  mv visit0000.png slidingSphareIB-cont.png
  cp ../../slidingSphareIB-iso.py .
  visit -cli -nowin -s slidingSphareIB-iso.py >& /dev/null
  mv visit0000.png slidingSphareIB-iso.png
endif
cd ../../..

echo '##################'
echo '### TENSION 2D ###'
echo '##################'
mkdir CheckVOF/Result_tmp/Tension2D
cp CheckVOF/main-surfaceTension-2D.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Tension2D
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpiexec -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
cat log.txt |grep velocity > vel.out
cp ../../vel.gnu .
gnuplot vel.gnu >& /dev/null
cp ../../tension2D.gth .
dogather

if ( $1 =~ "tec" ) then
  # tecplot
  cp ../../make_png.mcr .
  cp ../../tension2D.lay .
  tec360 tension2D.lay -b make_png.mcr >& /dev/null
  mv tmp.png tension2D.png
else
  # visit
  cp ../../tension2D.py .
  visit -cli -nowin -s tension2D.py >& /dev/null
  mv visit0000.png tension2D.png 
endif
cd ../../..

echo '##################'
echo '### TENSION 3D ###'
echo '##################'
mkdir CheckVOF/Result_tmp/Tension3D
cp CheckVOF/main-surfaceTension-3D.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Tension3D
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpiexec -np 16 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
cat log.txt |grep velocity > vel.out
cp ../../vel.gnu .
gnuplot vel.gnu >& /dev/null
cat log.txt |grep kappa_min_max > kappa.out
cp ../../kappa.gnu .
gnuplot kappa.gnu >& /dev/null
cp ../../tension3D.gth .
dogather
cp ../../tension3D2.gth . 
dogather

if ( $1 =~ "tec" ) then
  # tecplot
  cp ../../make_png.mcr .
  cp ../../tension3D.lay .
  tec360 tension3D.lay -b make_png.mcr >& /dev/null
  mv tmp.png tension3D.png
else
  # visit
  cp ../../tension3D.py .
  visit -cli -nowin -s tension3D.py >& /dev/null
  mv visit0000.png tension3D.png
endif
cd ../../..

echo '###############'
echo '### Enright ###'
echo '###############'
mkdir CheckVOF/Result_tmp/Enright
cp CheckVOF/main-enright.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Enright
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
if ( $2 =~ "sbatch" ) then
  cp ../../go-simple .
  sbatch go-simple
  echo '#\!/bin/tcsh' > post.sh
  echo 'cat log.txt |grep totalvol > vol.out' >> post.sh
  echo 'cp ../../vol.gnu .' >> post.sh
  echo 'gnuplot vol.gnu >& /dev/null' >> post.sh
  echo 'cp ../../enright.gth .' >> post.sh
  echo 'dogather' >> post.sh
  echo 'cp ../../enright.py .' >>post.sh
  echo 'visit -cli -nowin -s enright.py >& /dev/null' >>post.sh
  chmod 777 post.sh
else
  mpiexec -np 16 Boil > log.txt
  set sec1 = `date +%s`
  @ difft = $sec1 - $sec0
  echo $difft sec.
  # gnuplot
  cat log.txt |grep totalvol > vol.out
  cp ../../vol.gnu .
  gnuplot vol.gnu >& /dev/null
  cp ../../enright.gth .
  dogather

  if ( $1 =~ "tec" ) then
    # tecplot
    cp ../../make_png.mcr .
    cp ../../enright.lay .
    foreach i (*.plt)
      tec360 enright.lay -b make_png.mcr $i >& /dev/null
      mv tmp.png $i:r.png
    end
  else
    # visit
    cp ../../enright.py .
    visit -cli -nowin -s enright.py >& /dev/null
  endif
endif
cd ../../..

echo '###########################'
echo '### RISING BUBBLE IJNMF ###'
echo '###########################'
mkdir CheckVOF/Result_tmp/RisingBubbleIJNMF
cp CheckVOF/main-risingBubble-IJNMF.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/RisingBubbleIJNMF
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
if ( $2 =~ "sbatch" ) then
  cp ../../go-simple2 .
  sbatch go-simple2
  echo '#\!/bin/tcsh' > post.sh
  echo 'cat log.txt |grep totalvol > vol.out' >>post.sh
  echo 'cp ../../vol.gnu .' >>post.sh 
  echo 'gnuplot vol.gnu >& /dev/null' >>post.sh
  echo 'cat log.txt |grep front > front.out' >>post.sh
  echo 'cp ../../front-z.gnu .' >>post.sh
  echo 'gnuplot front-z.gnu >& /dev/null' >>post.sh
  echo 'cp ../../risingBubble-uvw.gth .' >>post.sh
  echo 'dogather' >>post.sh
  echo 'cp ../../risingBubble-xyz.gth .' >>post.sh
  echo 'dogather' >>post.sh
  echo 'cp ../../risingBubble-kappa.py .' >>post.sh
  echo 'visit -cli -nowin -s risingBubble-kappa.py >& /dev/null' >>post.sh
  echo 'mv visit0000.png risingBubble-kappa.png' >>post.sh
  echo 'cp ../../risingBubble.py .' >>post.sh
  echo 'visit -cli -nowin -s risingBubble.py >& /dev/null' >>post.sh
  echo 'mv visit0000.png risingBubble.png' >>post.sh
  echo 'cp ../../risingBubble-side.py .' >>post.sh
  echo 'visit -cli -nowin -s risingBubble-side.py >& /dev/null' >>post.sh
  chmod 777 post.sh
else
  mpiexec -np 16 Boil > log.txt
  set sec1 = `date +%s`
  @ difft = $sec1 - $sec0
  echo $difft sec.
  # gnuplot
  cat log.txt |grep totalvol > vol.out
  cp ../../vol.gnu .
  gnuplot vol.gnu >& /dev/null
  cat log.txt |grep front > front.out
  cp ../../front-z.gnu .
  gnuplot front-z.gnu >& /dev/null
  cp ../../risingBubble-uvw.gth .
  dogather
  cp ../../risingBubble-xyz.gth .
  dogather

  if ( $1 =~ "tec" ) then
    # tecplot
    cp ../../make_png.mcr .
    cp ../../risingBubble.lay .
    tec360 risingBubble.lay -b make_png.mcr >& /dev/null
    mv tmp.png risingBubble.png
    cp ../../risingBubble-kappa.lay .
    tec360 risingBubble-kappa.lay -b make_png.mcr >& /dev/null
    mv tmp.png risingBubble-kappa.png
    cp ../../risingBubble-side.lay .
    tec360 risingBubble-side.lay -b make_png.mcr >& /dev/null
    mv tmp.png risingBubble-side.png
  else
    # visit
    cp ../../risingBubble-kappa.py .
    visit -cli -nowin -s risingBubble-kappa.py >& /dev/null
    mv visit0000.png risingBubble-kappa.png
    cp ../../risingBubble.py .
    visit -cli -nowin -s risingBubble.py >& /dev/null
    mv visit0000.png risingBubble.png
    cp ../../risingBubble-side.py .
    visit -cli -nowin -s risingBubble-side.py >& /dev/null
  endif
endif
cd ../../..


exit
