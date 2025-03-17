D = 1.88e-9     # adjust manually, diffusion coefficient of species
dia = 0.009     # adjust manually
Ae = 3.1415 * dia**2

set term png
set gri
set xlabel "Time (s)"

set output "factor1.png"
set ylabel "Mass transfer factor (s^(-0.5))"
plot "masstransfer.out" u 2:(-$10/Ae/D**0.5) w l t "Method smdot"

set output "factor2.png"
plot "masstransfer.out" u 2:(-$12/Ae/D**0.5) w l t "Method flux"

