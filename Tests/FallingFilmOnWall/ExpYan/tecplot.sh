#! /usr/bin/tcsh
mkdir PLT
mkdir PNG-tpr
mkdir PNG-eps
mkdir PNG-mdot

foreach i (*.plt)
  tec360 $i tpr.lay -b ~/bin/make_png.mcr
  mv -f tmp.png PNG-tpr/$i:r.png

  tec360 $i eps.lay -b ~/bin/make_png.mcr
  mv -f tmp.png PNG-eps/$i:r.png

  tec360 $i mdot.lay -b ~/bin/make_png.mcr
  mv -f tmp.png PNG-mdot/$i:r.png

  mv $i PLT/.
end
