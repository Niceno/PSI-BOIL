set term png
set output "mdot.png"
set gri
set xlabel "Time (s)"
set ylabel "Mass transfer rate (kg/s)"
plot "mdot.out" u 2:4 w l t "Evaporation","" u 2:6 w l t "Condensation"
