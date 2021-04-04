      boil::print_line("######################");
      boil::print_line("#                    #");
      boil::print_line("# START FROM SCRATCH #");
      boil::print_line("#                    #");
      boil::print_line("######################");

      /* color */
      c = 1.0;
      conc.init();

      /* temperature */
      tpr = tout;
      for_vijk(tpr,i,j,k) {
        if(tpr.domain()->ibody().off(i,j,k)) {
          tpr[i][j][k] = twall;
        }
      }
      tpr.bnd_update();
      tpr.exchange_all();

      boil::plot->plot(uvw,c,tpr,press,
                       "uvw-c-tpr-press",
                       0);
