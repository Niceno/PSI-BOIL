cat log.txt |grep totalvol > vol.out
cat log.txt |grep ront_minmax > front.out
gnuplot vol.gnu
gnuplot front.gnu

