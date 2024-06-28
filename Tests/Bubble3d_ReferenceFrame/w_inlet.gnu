set term png
set gri
set xlabel "Time (s)"

set output "w_inlet.png"
set ylabel "Inlet velocity at top boundary (m/s)"
plot "w_inlet.out" u 2:4 w l

set output "accl_camera.png"
set ylabel "Acceleration of camera in Z (m/s^2)"
plot "w_inlet.out" u 2:6 w l

set output "pos_camera.png"
set ylabel "Position of camera in Z (m)"
plot "w_inlet.out" u 2:8 w l

set output "z_bubble.png"
set ylabel "Bubble top in Z (m)"
plot "w_inlet.out" u 2:10 w l
