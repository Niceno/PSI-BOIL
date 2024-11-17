set term png
set output "tpr-vfl.png"
set gri
set xlabel "X (m)"
set ylabel "Temperature (K)"
set y2tics
set y2label "Volume fraction of liquid"
plot "profiles-0.0012-000100.txt" u 1:6 w l t "Temperature","" u 1:5 w l axes x1y2 t "Volume fraction of liquid"
