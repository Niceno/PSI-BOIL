set term png
set output "dt.png"
set gri
set xlabel "Time step (-)"
set ylabel "Time increment (s)"
plot "dt.out" u 2 w l

