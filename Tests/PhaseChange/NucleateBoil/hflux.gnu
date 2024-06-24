set term png
set output "hflux.png"
set gri
set xlabel "Time (s)"
set ylabel "Heat flux (W)"
plot "hflux.out" u 2:3 w l t "total","" u 2:4 w l t "microlayer"

