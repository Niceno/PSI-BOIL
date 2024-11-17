cat log.txt | grep DT > dt.out
sed -i 's/#/ /g' dt.out

gnuplot dt.gnu

