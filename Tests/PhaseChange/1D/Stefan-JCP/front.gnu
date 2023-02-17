set term png
set output "front.png"
set xlabel "Time (s)"
set ylabel "Front position (m)"
set gri
plot 2*6.695e-2*sqrt(0.025*x/(0.597*2030)) ,"front.out" u 2:3 w l

