cat log.txt | grep DT > dt.out
cat log.txt | grep mdot_cutoff > mdot.out
sed -i 's/#/ /g' dt.out

gnuplot dt.gnu
gnuplot mdot.gnu


