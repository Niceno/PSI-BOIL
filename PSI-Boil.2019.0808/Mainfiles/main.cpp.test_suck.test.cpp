  real tval = tgamma(xpos);
  
  /* color function */
  c = 1.0;
  for_vijk(c,i,j,k) {
    if(c.xn(i+1)<xpos) {
      c[i][j][k]=0.0;
    } else if(c.xn(i)<xpos) {
      c[i][j][k]=1.0-(xpos-c.xn(i))/c.dxc(i);
    }
  }
  c.bnd_update();
  c.exchange_all();

  tpr = Tsat;
  /* temperature */
  for_vijk(tpr,i,j,k) {
    real xtmp=tpr.xc(i);
    if(xtmp>xpos) {
      tpr[i][j][k] = tsol(xtmp,tval); 
    }
  }
  tpr.bnd_update();
  tpr.exchange_all();

  conc.ancillary();
  conc.new_time_step();

  for_vijk(tpr,i,j,k) {
    real xtmp=tpr.xc(i);

    std::vector<real> values = {tpr[i][j][k],tpr[i-1][j][k],tpr[i+1][j][k]};
    std::vector<real> stencil = {0.,c.dxc(i),c.dxc(i)};
    if(cht.interface(Sign::neg(),Comp::i(),i,j,k)) {
      stencil[1] = cht.distance_int_x(Sign::neg(),i,j,k,values[1]);
    }
    real sd = (xtmp>xpos) ? secdiff(stencil,values) : 0.;

    boil::oout<<i<<" "<<c[i][j][k]<<" | "<<tpr[i][j][k]<<" "<<d1tdx(xtmp,tval)<<" "<<d2tdx(xtmp,tval)<<" | "
    <<cht.first_derivative(false,Comp::i(),i,j,k,AccuracyOrder::First(),false,Old::yes)<<" "
    <<cht.first_derivative(false,Comp::i(),i,j,k,AccuracyOrder::Second(),false,Old::yes)<<" "
    <<cht.first_derivative(false,Comp::i(),i,j,k,AccuracyOrder::Fourth(),false,Old::yes)<<" | "
    <<sd<<" "
    <<cht.second_derivative(false,Comp::i(),i,j,k,AccuracyOrder::Second(),false,Old::yes)<<" "
    <<cht.second_derivative(false,Comp::i(),i,j,k,AccuracyOrder::Fourth(),false,Old::yes)<<" | "
    <<boil::endl;
  }
  exit(0);

  /*---------------+
  |  phase change  |
  +---------------*/
  pc.update();
  real massflux_heat = pc.get_smdot();
  massflux_heat /= conc.topo->get_totarea();
  real dtdx_t2 = betasol*sqrt(alpv)/sqrt(time.current_time()+tval)*rhov;
  boil::oout<<"mflux= "<<tval<<" "
                       <<(xpos-x0)/DX<<" "
                       <<massflux_heat<<" "
                       <<dtdx_t2<<" "<<massflux_heat/dtdx_t2-1.
                       <<boil::endl;
