set term png
set output "w-vfl.png"
set gri
set xlabel "X (m)"
set ylabel "Velocity W (m/s)"
set y2tics
set y2label "Volume fraction of liquid"
plot "profiles-0.0012-000100.txt" u 1:4 w l t "W","" u 1:5 w l axes x1y2 t "Volume fraction of liquid"
