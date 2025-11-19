    /* test temperature */
    real tsol_max(-boil::unreal), tsol_min(boil::unreal);
    real tliq_max(-boil::unreal), tliq_min(boil::unreal);
    real tvap_max(-boil::unreal), tvap_min(boil::unreal);
    for_vijk(tpr,i,j,k) {
      real tval = tpr[i][j][k];
      if(tpr.domain()->ibody().on(i,j,k)) {
        if(cht.topo->above_interface(i,j,k)) {
          if(tval>tliq_max)
            tliq_max = tval;
          if(tval<tliq_min)
            tliq_min = tval;
        } else {
          if(tval>tvap_max)
            tvap_max = tval;
          if(tval<tvap_min)
            tvap_min = tval;
        }
      } else {
        if(tval>tsol_max)
          tsol_max = tval;
        if(tval<tsol_min)
          tsol_min = tval;
      }
    }
    boil::cart.max_real(&tsol_max);
    boil::cart.max_real(&tliq_max);
    boil::cart.max_real(&tvap_max);
    boil::cart.min_real(&tsol_min);
    boil::cart.min_real(&tliq_min);
    boil::cart.min_real(&tvap_min);
    boil::oout<<"tprextrema= "<<time.current_time()<<" "
              <<tsol_min<<" "<<tsol_max<<" "
              <<tliq_min<<" "<<tliq_max<<" "
              <<tvap_min<<" "<<tvap_max<<" "
              <<boil::endl;
    if(tliq_min<-0.1||tvap_min<-0.1||tsol_min<-0.1) {
      boil::oout<<"temperature instability. exiting."<<boil::endl;
      iint++;
      boil::plot->plot(uvw,c,tpr,mdot,mflx,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       iint);
      break;
    }
