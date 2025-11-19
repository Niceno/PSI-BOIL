set terminal png
set xlabel "Time (s)"
set ylabel "Front position (m)"
set gri
set output "front-z.png"
plot "front.out" u 2:4 w l

