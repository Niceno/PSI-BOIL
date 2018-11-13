set terminal png
set xlabel "Time (s)"
set ylabel "Max(|u|) and Max(|w|) (m/s)"
set gri
set output "vel.png"
plot "vel.out" u 2:3 w l,"" u 2:5 w l

