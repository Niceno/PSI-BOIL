rhov = 1.977   # adjust manually

set term png
set output "vol.png"

set gri
set xlabel "Time (s)"
set ylabel "Vapor volume (m3)"

set key left bottom

plot \
  "vol.out" using 2:4 w l t "Bubble volume", \
  "< paste vol.out mdot_x_dt_sum.out" using 2:($4-($6/rhov)) w l t "Bubble volume + Sum(mdot*dt)/rhov"

