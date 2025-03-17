dia = 0.009   # adjust manually
g = 9.8
T = sqrt(dia/g)

set term png
set gri

set output "sherwood.png"
set xlabel "Time (s)"
set ylabel "Sherwood nuber (-)"
plot "masstransfer.out" u 2:18 w l t "Method 1 (mdot)","" u 2:(-$20) w l t "Method 2 (flux at an elevation)"

set output "dia.png"
set ylabel "Diameter, de (m)"
plot "masstransfer.out" u 2:6 w l

set xlabel "Non-dimensional time, t'"

set output "re.png"
set ylabel "Re"
set y2label "Eo"
set y2tics
set yran [0:200]
plot "masstransfer.out" u ($2/T):16 axis x1y1 w l t "Re","" u ($2/T):14 axis x1y2 w l t "Eo"

set output "sh.png"
set ylabel "Sh"
set y2label "Eo"
set y2tics
set yran [0:10]
plot "masstransfer.out" u ($2/T):18 axis x1y1 w l t "Sh","" u ($2/T):14 axis x1y2 w l t "Eo"

set output "re-sh.png"
set xlabel "Re"
set ylabel "Sh"
unset y2tics
unset y2label
set autosc
set xran [0:]
plot "masstransfer.out" u 16:18 axis x1y1 w l

