    /* dynamic contact angle */
    int iglob, iproc;
    real clvel = conc.extract_cl_velocity_2d(Comp::i(),Comp::k(),Sign::neg(),
                                             &iglob,&iproc);
    real Cau = mul*std::max(boil::pico,clvel)/sig;

#if 0
    real theta3 = std::pow(cangle*boil::pi/180.,3.);
    real minval3 = std::pow(0.1*boil::pi/180.,3.);
    real cangnew = theta3-9.*Cau*log(0.5*dxmin/nanoscale);
    cangnew = 180./boil::pi*std::pow( std::max(minval3,cangnew), 1./3.);
#else
    real cangnew = std::min(Cae/Cau,thetahock)*180./boil::pi;
#endif
    conc.set_cangle(cangnew);
    boil::oout<<"dyncangle= "<<time.current_time()<<" "
              <<clvel<<" "<<cangnew<<" | "<<cangle<<" "<<iglob<<" "<<iproc<<boil::endl;
