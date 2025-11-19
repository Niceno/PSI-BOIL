      boil::print_line("######################");
      boil::print_line("#                    #");
      boil::print_line("# START FROM SCRATCH #");
      boil::print_line("#                    #");
      boil::print_line("######################");

      /* color */
      for(auto l : c.levels)
        *l = 1.0;

      conc_coarse.init();
      conc_fine.init();

      /* temperature */
      for(auto l : tpr.levels) {
        *l = tout;

        for_vijk((*l),i,j,k) {
          if(l->domain()->ibody().off(i,j,k)) {
            (*l)[i][j][k] = twall;
          }
        }
        l->bnd_update();
        l->exchange_all();
      }

      boil::plot->plot(uvw.coarse,c.coarse,tpr.coarse,press,
                       "uvw-c-tpr-press",
                       0);
