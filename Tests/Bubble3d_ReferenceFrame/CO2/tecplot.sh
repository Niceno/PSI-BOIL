#!/bin/tcsh
mkdir PNG-iso-eps
mkdir PNG-iso
mkdir PLT

foreach i (*.plt)
tec360 iso-eps.lay $i -b ~/bin/make_png.mcr
mv -f tmp.png ./PNG-iso-eps/$i:r.png
tec360 iso.lay $i -b ~/bin/make_png.mcr
mv -f tmp.png ./PNG-iso/$i:r.png
mv $i ./PLT/.
end

