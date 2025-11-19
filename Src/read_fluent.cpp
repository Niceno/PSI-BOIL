  /*----------------+
  | fluent-2DR3.txt |
  +----------------*/
  int nk_fluent;
  std::vector<real> z_fluent;
  std::vector<real> U_fluent;
  std::vector<real> tke_fluent;
  std::vector<real> omg_fluent;
  std::vector<real> len_fluent;
  std::vector<real> nut_fluent;
  std::fstream input_tke;
  input_tke.open("fluent-2DR3.txt", std::ios::in);
  if( !input_tke.fail() ) {
    input_tke >> nk_fluent;
    boil::oout<<"nk_fluent= "<<nk_fluent<<"\n";
    std::string line;
    input_tke >> line >> line >> line >> line >> line >> line >> line;
    for(int K=0; K<nk_fluent; K++) {
      real val1,val2,val3,val4,val5,val6,val7;
      input_tke >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
      z_fluent.push_back(val1);
      tke_fluent.push_back(val2);
      omg_fluent.push_back(val3);
      len_fluent.push_back(val5);
      U_fluent.push_back(val6);
      nut_fluent.push_back(val7);
    }
    //exit(0);
  } else {
    boil::oout<<"cannot find fluent-2DR3.txt"<<"\n";
    exit(0);
  }
  //for(int K=0; K<nk_fluent; K++) {
  //  boil::oout<<K<<" "<<z_fluent[K]<<" "<<tke_fluent[K]<<"\n";
  //}

#if 0
  // check how to access global
  boil::oout<<"nkg= "<<nkg<<" dom.gk "<<dom.gk()<<" dom.gik "<<dom.gik()<<"\n";
  // loop for cell
  for (int KC=0; KC<dom.gik(); KC++){
    boil::oout<<"cell: KC= "<<KC<<" dom.zc_global "<<dom.zc_global(KC+boil::BW)<<"\n";
  }
  // loop for node
  for (int KN=0; KN<=dom.gik(); KN++){
    boil::oout<<"node: KN= "<<KN<<" dom.zn_global "<<dom.zn_global(KN+boil::BW)<<"\n";
  }
#endif
  // store U, tke, turbulent length
  real U_global[dom.gik()+1], nut_global[dom.gik()], tke_global[dom.gik()],
       omg_global[dom.gik()], len_global[dom.gik()];
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if (ztmp<=0.0) {
      U_global[KC]=0.0;
      nut_global[KC]=0.0;
      tke_global[KC]=0.0;
      omg_global[KC]=0.0;
      len_global[KC]=0.0;
    } else {
      for(int KF=0; KF<nk_fluent; KF++) {
        if(ztmp<z_fluent[KF]) {
          real s1 = ztmp - z_fluent[KF-1];
          real s2 = z_fluent[KF] - ztmp;
          U_global[KC]   = s2/(s1+s2)*U_fluent[KF-1]  +s1/(s1+s2)*U_fluent[KF];
          nut_global[KC] = s2/(s1+s2)*nut_fluent[KF-1]+s1/(s1+s2)*nut_fluent[KF];
          tke_global[KC] = s2/(s1+s2)*tke_fluent[KF-1]+s1/(s1+s2)*tke_fluent[KF];
          omg_global[KC] = s2/(s1+s2)*omg_fluent[KF-1]+s1/(s1+s2)*omg_fluent[KF];
          len_global[KC] = s2/(s1+s2)*len_fluent[KF-1]+s1/(s1+s2)*len_fluent[KF];
          break;
        }
      }
    }
    //boil::oout<<"KC,Z,U,mut,tke,omg,len "<<KC<<" "<<dom.zc_global(KC+boil::BW)<<" "
    //          <<U_global[KC]<<" "<<nut_global[KC]<<" "
    //          <<tke_global[KC]<<" " <<omg_global[KC]<<" "<<len_global[KC]<<"\n";
  }
  // buffer cell
  U_global[dom.gik()]=U_global[dom.gik()-1];

  // calculate du/dz
  real dudz_global[dom.gik()];
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    real ztmp_minus = dom.zc_global(KC-1+boil::BW);
    if (ztmp<=0.0) {
      dudz_global[KC]=0.0;
    } else if (ztmp>0 && ztmp_minus<0) {
      dudz_global[KC]=(U_global[KC]-0.0)/(ztmp-0.0);
      //boil::oout<<"dudz at wall "<<U_global[KC]<<" "<<ztmp<<" "<<dudz_global[KC]<<" "<<KC<<"\n";
    } else {
      dudz_global[KC]=(U_global[KC+1]-U_global[KC-1])
                     /(dom.zc_global(KC+1+boil::BW)-dom.zc_global(KC-1+boil::BW));
      //if (KC==dom.gik()-1) {
      //  std::cout<<"ZC "<<dom.zc_global(KC+1+boil::BW)<<" "<<dom.zc_global(KC-1+boil::BW)<<"\n";
      //  exit(0);
      //}
    }
    //boil::oout<<"KC,ZC,U,dudz "<<KC<<" "<<dom.zc_global(KC+boil::BW)<<" "
    //          <<U_global[KC]<<" "<<dudz_global[KC]<<"\n";
  }

  // calculate Reynolds stress Rij(KC)
  real L11_global[dom.gik()], L22_global[dom.gik()], L31_global[dom.gik()], L33_global[dom.gik()];
  real sig_global[dom.gik()];
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if (ztmp<=0.0) {
      L11_global[KC]=0.0;
      L22_global[KC]=0.0;
      L31_global[KC]=0.0;
      L33_global[KC]=0.0;
    } else {
      real Rxx,Ryy,Rzz,Rxz,Rxz_max,dz;
      Rxx=2.0/3.0*tke_global[KC];
      Ryy=Rxx;
      Rzz=Rxx;
      Rxz=-nut_global[KC]*dudz_global[KC];
      Rxz_max = 0.98*sqrt(Rxx*Rzz);
      if (fabs(Rxz)>Rxz_max) Rxz=copysign(Rxz_max,Rxz);
      L11_global[KC]=sqrt(Rxx);
      L22_global[KC]=sqrt(Ryy);
      L31_global[KC]=Rxz/sqrt(Rxx);
      L33_global[KC]=sqrt(Rzz-Rxz*Rxz/Rxx);
      dz = dom.zn_global(KC+1+boil::BW)-dom.zn_global(KC+boil::BW);
      sig_global[KC]=max(dz,len_global[KC]);
      //boil::oout<<"KC,dz,len,sig "<<KC<<" "<<dz<<" "<<len_global[KC]<<" "<<sig_global[KC]<<"\n";
    }
  }

