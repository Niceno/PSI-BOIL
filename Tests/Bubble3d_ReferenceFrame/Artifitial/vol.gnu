rhov = 10.0    # adjust manually

set term png
set output "vol.png"

set gri
set xlabel "Time (s)"
set ylabel "Vapor volume (m3)"

set key left bottom

plot \
  "vol.out" using 2:4 w l t "Bubble volume"

