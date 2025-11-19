set terminal png
set xlabel "Time step"
set ylabel "Max and min kappa (1/m)"
set gri
set output "kappa.png"
plot "kappa.out" u 3 t "Min" w l,"" u 4 t "Max" w l

