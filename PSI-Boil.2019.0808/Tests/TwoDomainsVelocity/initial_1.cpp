
if( !restart ){

  boil::oout << "######################" << boil::endl;
  boil::oout << "# START FROM SCRATCH #" << boil::endl;
  boil::oout << "######################" << boil::endl;

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_avk(c_1,k){
    const real thickness = 0.005;
    const real wamp = 0.001;
    const real wave_length = LX / 4.0;
    const real wave_number = 2.0 * boil::pi / wave_length;
    if(c_1.zc(k)<(thickness-2.0*wamp)){
      for_avij(c_1,i,j){
        c_1[i][j][k]= 1.0;
      }
    } else if (c_1.zc(k)>(thickness+2.0*wamp)){
      for_avij(c_1,i,j){
        c_1[i][j][k]= 0.0;
      }
    } else {
      for_avij(c_1,i,j){
      int mm=8;
      real x0=dom_1.xn(i);
      real y0=dom_1.yn(j);
      real z0=dom_1.zn(k);
      real ddx=dom_1.dxc(i)/real(mm);
      real ddy=dom_1.dyc(j)/real(mm);
      real ddz=dom_1.dzc(k)/real(mm);
      int itmp=0;
      for (int ii=0; ii<mm; ii++){
        real xxc=x0+0.5*ddx+real(ii)*ddx;
        for (int jj=0; jj<mm; jj++){
          real wave = thickness + wamp * sin( wave_number * xxc);
          for (int kk=0; kk<mm; kk++){
            real zzc=z0+0.5*ddz+real(kk)*ddz;
            if ( zzc < wave){
              itmp=itmp+1;
            }
          }
        }
      }
      c_1[i][j][k]=real(itmp)/real(mm*mm*mm);
      }
    }
  }

  c_1.bnd_update();
  c_1.exchange_all();
  conc_1.init();
  update_step(c_1, step_1, sflag);

  for_vk(t_1,k) {
    if(c_1.zc(k)<0.000) {
      for_vij(c_1,i,j){
        t_1[i][j][k]=110.0;
      }
    } else {
      for_vij(c_1,i,j){
        t_1[i][j][k]=100.0;
      }
    }
  }
  t_1.bnd_update();
  t_1.exchange_all();

  boil::plot->plot(uvw_1,c_1,t_1,press_1,mu_t_1,
                  "d1-uvw-c-t-press-mu_t", 0);
}
