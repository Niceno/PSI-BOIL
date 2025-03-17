grep main:w_inlet log.txt > w_inlet.out
grep mass_transfer_rate:time log.txt > masstransfer.out
grep totalvol log.txt > vol.out
grep Dia,Eo,Re,Sh log.txt > dia_eo_re_sh.out
grep DT log.txt > dt.out
sed -i 's/# DT  : //' dt.out

paste masstransfer.out dt.out | awk '{
    product = $6 * $(NF)
    cum_sum += product
    print $2, cum_sum
}' > mdot_x_dt_sum.out


gnuplot w_inlet.gnu
gnuplot masstransfer.gnu
gnuplot factor.gnu
gnuplot dia_eo_re_sh.gnu
gnuplot vol.gnu

