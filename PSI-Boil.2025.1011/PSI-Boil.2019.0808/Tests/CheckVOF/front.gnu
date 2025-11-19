set terminal png
set xlabel "Time (s)"
set ylabel "Front position (m)"
set gri
set output "front.png"
plot "front.out" u 2:3 w l

