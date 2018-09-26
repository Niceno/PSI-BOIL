  
if( !restart ){

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vk(c_2,k) {
    if(c_2.zc(k)<0.000) {
      for_vij(c_2,i,j){
        c_2[i][j][k]=1.0;
        t_2[i][j][k]=110.0;
      }
    } else if(c_2.zc(k)<0.005) {
      for_vij(c_2,i,j){
        c_2[i][j][k]=1.0;
        t_2[i][j][k]=100.0;
      }
    } else {
      for_vij(c_2,i,j){
        c_2[i][j][k]=0.0;
        t_2[i][j][k]=100.0;
      }
    }
  }
  c_2.bnd_update();
  c_2.exchange_all();
  conc_2.init();
  update_step(c_2, step_2, sflag);

  t_2.exchange_all();

  boil::plot->plot(uvw_2,c_2,t_2,press_2,mu_t_2,mdot_2,
                  "d2-uvw-c-t-press-mu_t-mdot", 0);
}
