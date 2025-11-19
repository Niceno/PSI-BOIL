    /* monitoring point, every time step */
    boil::oout<<"m0.val= "<<time.current_time()<<" "<<time.dt()
              <<" U= "<<m0.value(uvw,Comp::u())<<" "
              <<m0.value(uvw,Comp::v())<<" "<<m0.value(uvw,Comp::w())<<"\n";
    /* monitoring every tint2 interval */
    if(int(time.current_time()/(tint2)) >= iint2 ) {
      iint2 = int(time.current_time() / tint2);
      boil::oout<<"main:Plotting:Rack\n";
      // r1
      r1.save(uvw,Comp::u(),"r1-U",iint2);
      r1.save(uvw,Comp::v(),"r1-V",iint2);
      r1.save(uvw,Comp::w(),"r1-W",iint2);
      r1.save(mu_t,"r1-mut",iint2);
      // r2
      r2.save(uvw,Comp::u(),"r2-U",iint2);
      r2.save(uvw,Comp::v(),"r2-V",iint2);
      r2.save(uvw,Comp::w(),"r2-W",iint2);
      r2.save(mu_t,"r2-mut",iint2);
      // r3
      r3.save(uvw,Comp::u(),"r3-U",iint2);
      r3.save(uvw,Comp::v(),"r3-V",iint2);
      r3.save(uvw,Comp::w(),"r3-W",iint2);
      r3.save(mu_t,"r3-mut",iint2);
      // r4
      r4.save(uvw,Comp::u(),"r4-U",iint2);
      r4.save(uvw,Comp::v(),"r4-V",iint2);
      r4.save(uvw,Comp::w(),"r4-W",iint2);
      r4.save(mu_t,"r4-mut",iint2);
      // r5
      r5.save(uvw,Comp::u(),"r5-U",iint2);
      r5.save(uvw,Comp::v(),"r5-V",iint2);
      r5.save(uvw,Comp::w(),"r5-W",iint2);
      r5.save(mu_t,"r5-mut",iint2);
      // increase counter
      // r0
      r0.save(uvw,Comp::u(),"r0-U",iint2);
      r0.save(uvw,Comp::v(),"r0-V",iint2);
      r0.save(uvw,Comp::w(),"r0-W",iint2);
      r0.save(mu_t,"r0-mut",iint2);
      // increase counter
      iint2 = int(time.current_time()/tint2) + 1;
    }
