cat log.txt |grep "twall=" > twall.out
cat log.txt |grep "plant_clr" > plant_clr.out
cat log.txt |grep "phasechange_update" > smdot.out
cat log.txt |grep "replant:sum_qsink_current_step" > qsink.out
cat log.txt |grep "enthalpyFD:hflux=" > hflux.out

gnuplot twall.gnu
gnuplot plant_clr.gnu
gnuplot qsink.gnu
gnuplot hflux.gnu


