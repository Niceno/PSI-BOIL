      boil::print_line("######################");
      boil::print_line("#                    #");
      boil::print_line("# START FROM SCRATCH #");
      boil::print_line("#                    #");
      boil::print_line("######################");

      /* color */
      real zcent, chord, V0;
      const real xcent = 0.0;
      boil::droplet_parameters_3D(180.-cangle,V0,const_cast<const real&>(radius),zcent,chord);

      boil::setup_circle_xz(conc.color(),radius,xcent,zcent);
      for_avijk(conc.color(),i,j,k)
        conc.color()[i][j][k] = 1. - conc.color()[i][j][k];

      conc.color().bnd_update();
      conc.color().exchange_all();
 
      conc.color_to_vf();
      conc.init();
      
      /* temperature */
      if(case_flag==1) {
        tpr = tout;

        for_vijk(tpr,i,j,k) {
          if(tpr.domain()->ibody().off(i,j,k)) {
            tpr[i][j][k] = twall;
          } 
        }
        tpr.bnd_update();
        tpr.exchange_all();
      } else {
        tpr = tout;

        interpolate_from_file(tpr);

        tpr.bnd_update();
        tpr.exchange_all();
      }

      boil::plot->plot(uvw,c,tpr,mdot,mflx,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       0);
