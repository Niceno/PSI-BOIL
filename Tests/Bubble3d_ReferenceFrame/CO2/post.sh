grep main:w_inlet log.txt > w_inlet.out
grep mass_transfer_rate:time log.txt > masstransfer.out
grep totalvol log.txt > vol.out

awk '{
    product = $10 * $4 
    cum_sum += product
    print $2, cum_sum
}' masstransfer.out > mdot_x_dt_sum.out


gnuplot w_inlet.gnu
gnuplot masstransfer.gnu
gnuplot factor.gnu
gnuplot vol.gnu

