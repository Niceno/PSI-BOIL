cat log.txt| grep Diameter= > diameter.out
cat log.txt| grep hflux: > hflux.out
cat log.txt |grep replant:site > t_nucl.out

gnuplot diameter.gnu
gnuplot t_nucl.gnu
gnuplot hflux.gnu

