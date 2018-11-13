set terminal png
set xlabel "Time (s)"
set ylabel "Volume (m^3)"
set gri
set output "vol.png"
plot "vol.out" u 2:3 w l

