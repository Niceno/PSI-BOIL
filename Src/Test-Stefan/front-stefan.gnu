set terminal png
set output "front-stefan.png"
set xlabel "Time (s)"
set ylabel "Front position (m)"
set gri
plot "front.out" u 2:3 w l

