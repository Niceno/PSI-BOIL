set term png
set output "eps-vfl.png"
set gri
set xlabel "X (m)"
set ylabel "Volume fraction of air in gas, eps"
set y2tics
set y2label "Volume fraction of liquid"
plot "profiles-0.0012-000100.txt" u 1:7 w l t "eps","" u 1:5 w l axes x1y2 t "Volume fraction of liquid"
