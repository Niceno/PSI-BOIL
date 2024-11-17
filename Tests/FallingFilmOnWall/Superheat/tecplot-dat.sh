#! /usr/bin/tcsh
mkdir PLT
mkdir PNG-tpr
mkdir PNG-eps

foreach i (*.dat)
  preplot $i

  tec360 $i:r.plt tpr.lay -b ~/bin/make_png.mcr
  mv -f tmp.png PNG-tpr/$i:r.png

  tec360 $i:r.plt eps.lay -b ~/bin/make_png.mcr
  mv -f tmp.png PNG-eps/$i:r.png

  mv $i:r.plt PLT/.
  rm $i
end
