set term png
set output "diameter.png"
set gri
set xlabel "Time (s)"
set ylabel "Diameter (m)"
plot "diameter.out" u 2:3 w l

