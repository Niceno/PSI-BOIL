set terminal png
set output "t_nucl.png"
set gri
set xlabel "Time (s)" 
set ylabel "Temperature at nucleation site (deg.)" 
plot "t_nucl.out" u 2:3 w l
