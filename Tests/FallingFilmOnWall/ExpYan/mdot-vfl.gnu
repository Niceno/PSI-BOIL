set term png
set output "mdot-vfl.png"
set gri
set xlabel "X (m)"
set ylabel "Mass transfer rate [kg/sm^3], mdot"
set y2tics
set y2label "Volume fraction of liquid"
plot "profiles-0.0012-000100.txt" u 1:8 w l t "mdot","" u 1:5 w l axes x1y2 t "Volume fraction of liquid"
