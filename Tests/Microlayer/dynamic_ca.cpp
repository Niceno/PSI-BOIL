    /* dynamic contact angle */
    real theta3 = std::pow(cangle*boil::pi/180.,3.);
    real minval3 = std::pow(0.1*boil::pi/180.,3.);
    real nanoscale = 10e-9;
    int iglob, iproc;
    real clvel = conc.extract_cl_velocity_2d(Comp::i(),Comp::k(),Sign::neg(),
                                             &iglob,&iproc);
    real cangnew = theta3-9.*mul*std::max(0.0,clvel)/sig*log(0.5*dxmin/nanoscale);
    cangnew = 180./boil::pi*std::pow( std::max(minval3,cangnew), 1./3.);
    conc.set_cangle(cangnew);
    boil::oout<<"dyncangle= "<<time.current_time()<<" "
              <<clvel<<" "<<cangnew<<" | "<<cangle<<" "<<iglob<<" "<<iproc<<boil::endl;
