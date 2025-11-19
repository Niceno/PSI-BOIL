set terminal png
set output "front-sucking.png"
set xlabel "Time (s)"
set ylabel "Front position (m)"
set gri
xint = 2e-4
t0 = 0.1
plot "front-sucking.out" u ($2+t0):($3-xint) w l

