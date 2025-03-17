dia = 2.176e-4
g = 9.8
T = sqrt(dia/g)

set term png
set grid
set xlabel "Non-dimensional time, t'"

set output "dia.png"
set ylabel "Diameter (m)
set y2label "Eo"
set y2tics
plot "dia_eo_re_sh.out" u ($2/T):3 axis x1y1 w l t "diameter","" u ($2/T):4 axis x1y2 w l t "Eo"

set output "re.png"
set ylabel "Re"
set y2label "Eo"
set y2tics
plot "dia_eo_re_sh.out" u ($2/T):5 axis x1y1 w l t "Re","" u ($2/T):4 axis x1y2 w l t "Eo"

set output "sh.png"
set ylabel "Sh"
set y2label "Eo"
set y2tics
set xran [0:20]
set yran [0:10]
set y2ran [0:2.5]
plot "dia_eo_re_sh.out" u ($2/T):6 axis x1y1 w l t "Sh","" u ($2/T):4 axis x1y2 w l t "Eo"

