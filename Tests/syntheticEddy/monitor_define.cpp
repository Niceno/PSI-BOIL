  /*--------------------+
  |  monitoring points  |
  +--------------------*/
  // find I
  real Xr1=0.001, Xr2=0.003, Xr3=0.005, Xr4=0.007, Xr5=0.009, Xr0=dx*0.5;
  int Ir1,Ir2,Ir3,Ir4,Ir5,Ir0;
  for(int I=boil::BW; I<=NX1+boil::BW-1; I++) {
    real xtmp = dom.xc_global(I);
    //boil::oout<<"I "<<I<<" xtmp "<<xtmp<<" "<<dx<<"\n";
    if(Xr1<=xtmp && xtmp<Xr1+dx) Ir1=I-boil::BW;
    if(Xr2<=xtmp && xtmp<Xr2+dx) Ir2=I-boil::BW;
    if(Xr3<=xtmp && xtmp<Xr3+dx) Ir3=I-boil::BW;
    if(Xr4<=xtmp && xtmp<Xr4+dx) Ir4=I-boil::BW;
    if(Xr5<=xtmp && xtmp<Xr5+dx) Ir5=I-boil::BW;
    //if(Xr0<=xtmp && xtmp<Xr0+dx) Ir0=I-boil::BW;
  }
  Ir0=1;
  boil::oout<<"main:Ir1,Ir2,Ir3,Ir4,Ir5,Ir0="<<Ir1<<" "<<Ir2<<" "<<Ir3<<" "<<Ir4<<" "<<Ir5<<" "<<Ir0<<"\n";
  Location m0("m0", dom, Ir1+1, NY/2, 4);  // 4 cell above wall
  // x~0.001
  Rack r1("r1", dom, Ir1+1, Range<int>(1,NY), Range<int>(1,NZ));
  for_vijk(p,i,j,k){p[i][j][k]=p.dV(i,j,k);}
  r1.save_grid(p,"r2-scalar");
  r1.save(p,"r1-volume",0);
  r1.save_grid(uvw,Comp::u(),"r1-vector-U");
  r1.save_grid(uvw,Comp::v(),"r1-vector-V");
  r1.save_grid(uvw,Comp::w(),"r1-vector-W");

  // x~0.003
  Rack r2("r2", dom, Ir2, Range<int>(1,NY), Range<int>(1,NZ));
  r2.save_grid(p,"r2-scalar");
  r2.save(p,"r2-volume",0);
  r2.save_grid(uvw,Comp::u(),"r2-vector-U");
  r2.save_grid(uvw,Comp::v(),"r2-vector-V");
  r2.save_grid(uvw,Comp::w(),"r2-vector-W");

  // x~0.005
  Rack r3("r3", dom, Ir3, Range<int>(1,NY), Range<int>(1,NZ));
  r3.save_grid(p,"r3-scalar");
  r3.save(p,"r3-volume",0);
  r3.save_grid(uvw,Comp::u(),"r3-vector-U");
  r3.save_grid(uvw,Comp::v(),"r3-vector-V");
  r3.save_grid(uvw,Comp::w(),"r3-vector-W");

  // x~0.007
  Rack r4("r4", dom, Ir4, Range<int>(1,NY), Range<int>(1,NZ));
  r4.save_grid(p,"r4-scalar");
  r4.save(p,"r4-volume",0);
  r4.save_grid(uvw,Comp::u(),"r4-vector-U");
  r4.save_grid(uvw,Comp::v(),"r4-vector-V");
  r4.save_grid(uvw,Comp::w(),"r4-vector-W");

  // x~0.009
  Rack r5("r5", dom, Ir5, Range<int>(1,NY), Range<int>(1,NZ));
  r5.save_grid(p,"r5-scalar");
  r5.save(p,"r5-volume",0);
  p=0.0;
  r5.save_grid(uvw,Comp::u(),"r5-vector-U");
  r5.save_grid(uvw,Comp::v(),"r5-vector-V");
  r5.save_grid(uvw,Comp::w(),"r5-vector-W");

  // x~0.00
  Rack r0("r0", dom, Ir0, Range<int>(1,NY), Range<int>(1,NZ));
  r0.save_grid(p,"r0-scalar");
  r0.save(p,"r0-volume",0);
  p=0.0;
  r0.save_grid(uvw,Comp::u(),"r0-vector-U");
  r0.save_grid(uvw,Comp::v(),"r0-vector-V");
  r0.save_grid(uvw,Comp::w(),"r0-vector-W");

