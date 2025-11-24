set gri
set xlabel "Time (s)"
set ylabel "Radius (m)"
beta = 4.06022
Cp = 4215.9*958.4
lambda = 0.679
r0 = 50e-6
t0 = r0**2.0*Cp/(4.0*lambda*beta**2.0)

set term "png"
set output "radius.png"
plot "front.out" u ($2+t0):4 w l,2.0*beta*sqrt(lambda/Cp*x) w l
