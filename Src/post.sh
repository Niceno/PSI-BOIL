cat log.txt | grep "replant:time" > replant.out
cat log.txt | grep "totalvol" > vol.out
gnuplot replant.gnu
gnuplot vol.gnu

