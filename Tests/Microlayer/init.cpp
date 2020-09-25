  int ts;
  /* load variables */
  std::vector<Scalar*> load_scalars = { &press, &c.coarse };
  if(case_flag!=0)
    load_scalars.push_back(&tpr.fine);
  else
    load_scalars.push_back(&tpr.coarse);
  std::vector<std::string> load_scalar_names = { "press", "c", "tpr" };

  std::vector<Vector*> load_vectors = { &uvw.coarse };
  std::vector<std::string> load_vector_names = { "uvw" };

  /* solid */
  csol.fine   = 0.0;
  csol.coarse = 0.0;
  for(auto l : csol.levels) {
    for_vijk((*l),i,j,k) {
      if(l->zn(k+1)<LZheat) {
        (*l)[i][j][k] = 1.0;
      } else if(l->zn(k)<= -LZheat) {
        (*l)[i][j][k] = (LZheat+l->zn(k+1))/l->dzc(k);
      }
    }
  }

  if(boil::load_backup("time.txt",ts,time,
                       load_scalars, load_scalar_names,
                       load_vectors, load_vector_names)) {
    conc_coarse.init();
    boil::prolongate_color_XZ(conc_coarse,conc_fine);
    conc_fine.color_to_vf();
    conc_fine.init();
  } else {

    /* start from single phase scratch */
    if(case_flag==0) {
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

    /* create a bubble from initial time step */
    } else if(case_flag==1) {
      boil::print_line("######################");
      boil::print_line("#                    #");
      boil::print_line("# START FROM SCRATCH #");
      boil::print_line("#                    #");
      boil::print_line("######################");

      /* color */
      real zcent, chord, V0;
      const real xcent = 0.0;
      boil::droplet_parameters_3D(180.-cangle,V0,const_cast<const real&>(radius),zcent,chord);

      boil::setup_circle_xz(conc_fine.color(),radius,xcent,zcent);
      for_avijk(conc_fine.color(),i,j,k)
        conc_fine.color()[i][j][k] = 1. - conc_fine.color()[i][j][k];

      conc_fine.color().bnd_update();
      conc_fine.color().exchange_all();
 
      conc_fine.color_to_vf();
      conc_fine.init();
      
      c.restrict_volume_XZ();
      conc_coarse.init();

      /* temperature */
      for_coarsefine(l) {
        tpr[l] = tout;

        for_vijk(tpr[l],i,j,k) {
          if(tpr[l].domain()->ibody().off(i,j,k)) {
            tpr[l][i][j][k] = twall;
          } else if(c[l][i][j][k]<0.5) { /* almost good */
            tpr[l][i][j][k] = tsat0;
          } else if(tpr[l].zc(k)<=ztconst) {
            tpr[l][i][j][k] = twall + (tout-twall)/ztconst * tpr[l].zc(k);
          }
        }
        tpr[l].bnd_update();
        tpr[l].exchange_all();
      }

      boil::plot->plot(uvw.coarse,c.coarse,press,
                       "uvwc-cc-press",
                       0);

      boil::plot->plot(uvw.fine,conc_fine.color(),tpr.fine,mdot,mflx,
                       "uvwf-clrf-tprf-mdot-mflx",
                       0);

    } else {
      OMS(Underdevelopment. Exiting.);
      exit(0);
    }

  }

  /* set iint */
  int iint = time.current_time() / t_per_plot;
  boil::oout<<"iint= "<<iint<<"\n";
