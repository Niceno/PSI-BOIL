#if 1
  // calculate len_ave
  real len_ave=0.0;
  real sum_dz=0.0;
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if(ztmp>0.0){
      //real dz = dom.dzc_global(KC+boil::BW);
      real dz = dom.zn_global(KC+boil::BW+1)-dom.zn_global(KC+boil::BW);
      len_ave += dz*len_global[KC];
      sum_dz += dz;
    }
  }
  len_ave = len_ave/sum_dz;  // average sigma
  boil::aout<<"main:len_ave= "<<len_ave<<"\n";
  for (int KC=0; KC<dom.gik(); KC++) {
    len_global[KC]=len_ave;
  }
#endif

  // calculate TKE and turbulent length at cell center
  real tke_local[p.ek()+1], len_local[p.ek()+1], u_local[p.ek()+1];
  for_vk(p,k) {
    //boil::oout<<"k= "<<k<<"\n";
    real ztmp = p.zc(k);
    if (ztmp<=0.0) {
      tke_local[k]=0.0;
      len_local[k]=0.0;
      u_local[k]=0.0;
    } else {
      for(int K=0; K<nk_fluent; K++) {
        if(ztmp<z_fluent[K]) {
          real s1 = ztmp - z_fluent[K-1];
          real s2 = z_fluent[K] - ztmp;
          real len_tmp;
          tke_local[k] = s2/(s1+s2)*tke_fluent[K-1]+s1/(s1+s2)*tke_fluent[K];
          len_local[k] = s2/(s1+s2)*len_fluent[K-1]+s1/(s1+s2)*len_fluent[K];
          u_local[k] = s2/(s1+s2)*U_fluent[K-1]+s1/(s1+s2)*U_fluent[K];
          //continue;
          break;
        }
      }
    }
  }
#if 1
  for_vk(p,k) {
    boil::oout<<k<<" z "<<p.zc(k)<<" tke "<<tke_local[k]<<" len "<<len_local[k]<<" U "<<u_local[k]<<"\n";
  }
#endif
  // calculate TKE and turbulent length at w-cell center
  real tke_local_w[p.ek()+1], len_local_w[p.ek()+1], u_local_w[p.ek()+1];
  for_vmk(uvw,Comp::w(),k) {
    //boil::oout<<"k= "<<k<<"\n";
    real ztmp = uvw.zc(Comp::w(),k);
    if (ztmp<=0.0) {
      tke_local_w[k]=0.0;
      len_local_w[k]=0.0;
      u_local_w[k]=0.0;
    } else {
      for(int K=0; K<nk_fluent; K++) {
        if(ztmp<z_fluent[K]) {
          real s1 = ztmp - z_fluent[K-1];
          real s2 = z_fluent[K] - ztmp;
          tke_local_w[k] = s2/(s1+s2)*tke_fluent[K-1]+s1/(s1+s2)*tke_fluent[K];
          len_local_w[k] = s2/(s1+s2)*len_fluent[K-1]+s1/(s1+s2)*len_fluent[K];
          u_local_w[k]   = s2/(s1+s2)*U_fluent[K-1]+s1/(s1+s2)*U_fluent[K];
          //continue;
          break;
        }
      }
    }
  }
#if 1
  for_vk(p,k) {
    boil::oout<<k<<" "<<uvw.zc(Comp::w(),k)<<" "<<tke_local_w[k]<<" "<<len_local_w[k]<<" "<<u_local_w[k]<<"\n";
  }
#endif

