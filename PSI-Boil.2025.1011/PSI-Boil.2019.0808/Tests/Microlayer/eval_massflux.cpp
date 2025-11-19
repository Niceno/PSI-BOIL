    real massflow_heat = pc.get_smdot();
    real massflux_heat = massflow_heat/conc.topo->get_totarea();
    real massflux_inert = rhov*sqrt(boil::pi/7.*rhov*latent*deltat_nucl
                                    /rhol/tsat0_K);

    /* transition zone */
    real xmaxbub = conc.topo->get_xmaxft();
    real xmaxml = output_row(ml_range_z.last(),c,true);
    real xtrml = 1.5*xmaxml - 0.5*xmaxbub;
    boil::oout<<"MLranges= "<<time.current_time()<<" "<<xtrml
                            <<" "<<xmaxml<<" "<<xmaxbub<<boil::endl;
    Range<real> ml_range_x(xtrml,xmaxml);

    real massflow_cap(0.0), massflow_ml(0.0);
    real are_cap(0.0), are_ml(0.0);
    for_vijk(c,i,j,k) {
      if(conc.topo->interface(i,j,k)) {
        real a = conc.topo->get_area(i,j,k);
        real a_cap = cap_frac(c.xc(i),
                              c.zc(k),
                              ml_range_x, ml_range_z,
                              a);
        real a_ml = a-a_cap;  

        massflow_cap += mflx[i][j][k]*a_cap;
        massflow_ml += mflx[i][j][k]*a_ml;
        are_cap += a_cap;
        are_ml += a_ml;
      }
    }
    boil::cart.sum_real(&massflow_cap);
    boil::cart.sum_real(&massflow_ml);
    boil::cart.sum_real(&are_cap);
    boil::cart.sum_real(&are_ml);
    
    real massflux_cap(0.0), massflux_ml(0.0);
    if(are_cap>0.0)
      massflux_cap = massflow_cap/are_cap;
    if(are_ml>0.0)
      massflux_ml = massflow_ml/are_ml;

#if 0
    for_vijk(c,i,j,k) {
      if(conc.topo->interface(i,j,k)) {
        mflx[i][j][k] = cap_val(c.xc(i),
                                c.zc(k),
                                ml_range_x, ml_range_z,
                                mflx[i][j][k],
                                massflux_cap);
      }
    }
    mflx.bnd_update();
    mflx.exchange();
    pc.finalize();
#elif 0
    mflx = massflux_heat;
    mflx.bnd_update();
    mflx.exchange();
    pc.finalize();
#endif

    boil::oout<<"mflux= "<<time.current_time()<<" "
                         <<massflux_heat<<" "<<massflux_inert<<" "
                         <<massflux_inert/massflux_heat<<" | "
                         <<massflow_heat<<" "<<massflow_cap<<" "<<massflow_ml<<" | "
                         <<massflux_cap<<" "<<massflux_ml
                         <<boil::endl;
