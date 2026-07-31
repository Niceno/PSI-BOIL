#include "vof.h"

#if 1
/******************************************************************************/
void VOF::front_minmax() {
  front_minmax(Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa) );
}

/******************************************************************************/
void VOF::front_minmax(Range<real> xr,
                       Range<real> yr,
                       Range<real> zr ) {
/***************************************************************************//**
*  \brief Detect maximum and minimum of free surface position.
*           results: xminft,xmaxft,yminft,ymaxft,zminft,zmaxft
*******************************************************************************/

   xminft=boil::exa; xmaxft=-boil::exa;
   yminft=boil::exa; ymaxft=-boil::exa;
   zminft=boil::exa; zmaxft=-boil::exa;
   //bool frontExist=false;

   real xfront, yfront, zfront;
   real phim, phip;

   /* i-direction */
   for(int i=si(); i<=ei()+1; i++) {
     if (phi.xc(i  )<xr.first()) continue;
     if (phi.xc(i+1)>xr.last() ) continue;
     for(int j=sj(); j<=ej()  ; j++) {
       if (phi.yc(j  )<yr.first()) continue;
       if (phi.yc(j  )>yr.last() ) continue;
       for(int k=sk(); k<=ek()  ; k++) {
         if (phi.zc(k  )<zr.first()) continue;
         if (phi.zc(k  )>zr.last() ) continue;
         if ( dom->ibody().off(i,j,k) || dom->ibody().off(i+1,j,k)) continue;
         phim=phi[i-1][j][k];
         phip=phi[i  ][j][k];
         if((phim-phisurf)*(phip-phisurf)<=0.0) {
           //frontExist=true;
           xfront=frontPosition(i,j,k,Comp::i());
           yfront=phi.yc(j);
           zfront=phi.zc(k);
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) yminft=yfront;
           if(yfront>ymaxft) ymaxft=yfront;
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
         }
       }  /* k */
     }  /* j */
   } /* i */

   /* j-direction */
   for(int i=si(); i<=ei()  ; i++) {
     if (phi.xc(i  )<xr.first()) continue;
     if (phi.xc(i  )>xr.last() ) continue;
     for(int j=sj(); j<=ej()+1; j++) {
       if (phi.yc(j  )<yr.first()) continue;
       if (phi.yc(j+1)>yr.last() ) continue;
       for(int k=sk(); k<=ek()  ; k++) {
         if (phi.zc(k  )<zr.first()) continue;
         if (phi.zc(k  )>zr.last() ) continue;
         if ( dom->ibody().off(i,j,k) || dom->ibody().off(i,j+1,k)) continue;
         phim=phi[i][j-1][k];
         phip=phi[i][j  ][k];
         if((phim-phisurf)*(phip-phisurf)<=0.0) {
           //frontExist=true;
           xfront=phi.xc(i);
           yfront=frontPosition(i,j,k,Comp::j());
           zfront=phi.zc(k);
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) yminft=yfront;
           if(yfront>ymaxft) ymaxft=yfront;
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
         } 
       } /* k */
     } /* j */
   } /* i */

   /* k-direction */
   for(int i=si(); i<=ei()  ; i++) {
     if (phi.xc(i  )<xr.first()) continue;
     if (phi.xc(i  )>xr.last() ) continue;
     for(int j=sj(); j<=ej()  ; j++) {
       if (phi.yc(j  )<yr.first()) continue;
       if (phi.yc(j  )>yr.last() ) continue;
       for(int k=sk(); k<=ek()+1; k++) {
         if (phi.zc(k  )<zr.first()) continue;
         if (phi.zc(k+1)>zr.last() ) continue;
         if ( dom->ibody().off(i,j,k) || dom->ibody().off(i,j,k+1)) continue;
         phim=phi[i][j][k-1];
         phip=phi[i][j][k  ];
         if((phim-phisurf)*(phip-phisurf)<=0.0) {
           //frontExist=true;
           xfront=phi.xc(i);
           yfront=phi.yc(j);
           zfront=frontPosition(i,j,k,Comp::k());
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) yminft=yfront;
           if(yfront>ymaxft) ymaxft=yfront;
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
         }
       } /* k */
     } /* j */
   } /* i */

   boil::cart.min_real(&xminft);
   boil::cart.max_real(&xmaxft);
   boil::cart.min_real(&yminft);
   boil::cart.max_real(&ymaxft);
   boil::cart.min_real(&zminft);
   boil::cart.max_real(&zmaxft);

   std::cout.setf(std::ios_base::scientific);
   boil::oout << "front_minmax:time,xmax= " << time->current_time() 
              << " " << xmaxft << " " << zmaxft << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);
}

real VOF::frontPosition(const int i, const int j, const int k,
                        const Comp m){
   real xyzfront;
   if (m == Comp::i()) { 
     xyzfront = fs[m][i][j][k];
     if(!boil::realistic(xyzfront))
       xyzfront = phi.xn(i);
   } else if (m == Comp::j()) {
     xyzfront = fs[m][i][j][k];
     if(!boil::realistic(xyzfront))
       xyzfront = phi.yn(j);
   } else {
     xyzfront = fs[m][i][j][k];
     if(!boil::realistic(xyzfront))
       xyzfront = phi.zn(k);
   }

   return xyzfront;
}
#endif
#if 0
real frontPositionIn(real xyz1, real xyz2, real phi1, real phi2);

/******************************************************************************/
void VOF::front_minmax() {
  front_minmax(Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa) );
}
/******************************************************************************/
void VOF::front_minmax(Range<real> xr,
                           Range<real> yr,
                           Range<real> zr ) {
/***************************************************************************//**
*  \brief Detect maximu and minimum of free surface position.
*           results: xminft,xmaxft,yminft,ymaxft,zminft,zmaxft
*******************************************************************************/

  xminft=boil::exa; xmaxft=-boil::exa;
  yminft=boil::exa; ymaxft=-boil::exa;
  zminft=boil::exa; zmaxft=-boil::exa;

  real phi1, phi2, xyz1, xyz2, xfront, yfront, zfront;

  real ymin_tmp[3],ymax_tmp[3];
  int ymax_ijk[3];

  for (int i=phi.si()-1; i<=phi.ei()  ; i++) {
    //if (phi.xc(i  )<xr.first()) continue;
    //if (phi.xc(i+1)>xr.last() ) continue;
    for (int j=phi.sj()  ; j<=phi.ej()  ; j++) {
      //if (phi.yc(j  )<yr.first()) continue;
      //if (phi.yc(j  )>yr.last() ) continue;
      for (int k=phi.sk()  ; k<=phi.ek()  ; k++) {
        //if (phi.zc(k  )<zr.first()) continue;
        //if (phi.zc(k  )>zr.last() ) continue;
        /* i-direction */
        phi1=phi[i  ][j][k];
        phi2=phi[i+1][j][k];
#ifdef IB
        if ( dom->ibody().off(i,j,k) || dom->ibody().off(i+1,j,k)) continue;
#endif
        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           xyz1=phi.xc(i);
           xyz2=phi.xc(i+1);
           xfront=frontPositionIn(xyz1,xyz2,phi1,phi2);
           yfront=phi.yc(j);
           zfront=phi.zc(k);
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           //if(yfront<yminft) yminft=yfront;
           //if(yfront>ymaxft) ymaxft=yfront;
           //if(zfront<zminft) zminft=zfront;
           //if(zfront>zmaxft) zmaxft=zfront;
        }
      }
    }
  }
  for (int i=phi.si()  ; i<=phi.ei()  ; i++) {
    //if (phi.xc(i  )<xr.first()) continue;
    //if (phi.xc(i  )>xr.last() ) continue;
    for (int j=phi.sj()-1; j<=phi.ej()  ; j++) {
      //if (phi.yc(j  )<yr.first()) continue;
      //if (phi.yc(j+1)>yr.last() ) continue;
      for (int k=phi.sk()  ; k<=phi.ek()  ; k++) {
        //if (phi.zc(k  )<zr.first()) continue;
        //if (phi.zc(k  )>zr.last() ) continue;
        /* j-direction */
        phi1=phi[i][j  ][k];
        phi2=phi[i][j+1][k];
#ifdef IB
        if ( dom->ibody().off(i,j,k) || dom->ibody().off(i,j+1,k)) continue;
#endif
        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           xyz1=phi.yc(j);
           xyz2=phi.yc(j+1);
           xfront=phi.xc(i);
           yfront=frontPositionIn(xyz1,xyz2,phi1,phi2);
           zfront=phi.zc(k);
           //if(xfront<xminft) xminft=xfront;
           //if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) {
             yminft=yfront;
             ymin_tmp[0]=xfront;
             ymin_tmp[1]=yfront;
             ymin_tmp[2]=zfront;
           }
           if(yfront>ymaxft) {
             ymaxft=yfront;
             ymax_tmp[0]=xfront;
             ymax_tmp[1]=yfront;
             ymax_tmp[2]=zfront;
             ymax_ijk[0]=i;
             ymax_ijk[1]=j;
             ymax_ijk[2]=k;
           }
           //if(zfront<zminft) zminft=zfront;
           //if(zfront>zmaxft) zmaxft=zfront;
        }
      }
    }
  }
  for (int i=phi.si()  ; i<=phi.ei()  ; i++) {
    //if (phi.xc(i  )<xr.first()) continue;
    //if (phi.xc(i  )>xr.last() ) continue;
    for (int j=phi.sj()  ; j<=phi.ej()  ; j++) {
      //if (phi.yc(j  )<yr.first()) continue;
      //if (phi.yc(j  )>yr.last() ) continue;
      for (int k=phi.sk()-1; k<=phi.ek()  ; k++) {
        //if (phi.zc(k  )<zr.first()) continue;
        //if (phi.zc(k+1)>zr.last() ) continue;
        /* k-direction */
        phi1=phi[i][j][k  ];
        phi2=phi[i][j][k+1];
#ifdef IB
        if ( dom->ibody().off(i,j,k) || dom->ibody().off(i,j,k+1)) continue;
#endif
        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           xyz1=phi.zc(k);
           xyz2=phi.zc(k+1);
           xfront=phi.xc(i);
           yfront=phi.yc(j);
           zfront=frontPositionIn(xyz1,xyz2,phi1,phi2);
           //if(xfront<xminft) xminft=xfront;
           //if(xfront>xmaxft) xmaxft=xfront;
           //if(yfront<yminft) yminft=yfront;
           //if(yfront>ymaxft) ymaxft=yfront;
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
         }
       }
     }
   }
   boil::cart.min_real(&xminft);
   boil::cart.max_real(&xmaxft);
   boil::cart.min_real(&yminft);
   boil::cart.max_real(&ymaxft);
   boil::cart.min_real(&zminft);
   boil::cart.max_real(&zmaxft);

   std::cout.setf(std::ios_base::scientific);
   boil::oout << "vof_front_minmax:time,xmax= " << time->current_time()
              << " " << xmaxft << " " << zmaxft << boil::endl;
   boil::oout << "ymax_tmp= "<<ymaxft<<" "<<ymax_tmp[0]<<" "<<ymax_tmp[1]<<" "<<ymax_tmp[2]<<"\n";
   boil::oout << "ymax_ijk= "<<ymax_ijk[0]<<" "<<ymax_ijk[1]<<" "<<ymax_ijk[2]<<"\n";
   boil::oout << "VOF(56,58,24):i "<<phi[55][58][24]<<" "<<phi[56][58][24]<<" "<<phi[57][58][24]<<"\n";
   boil::oout << "VOF(56,58,24):j "<<phi[56][57][24]<<" "<<phi[56][58][24]<<" "<<phi[56][59][24]<<"\n";
   boil::oout << "VOF(56,58,24):k "<<phi[56][58][23]<<" "<<phi[56][58][24]<<" "<<phi[56][58][25]<<"\n";
   boil::oout << "VOF(37,56,35)= "<<phi[37][55][35]<<" "<<phi[37][56][35]<<" "<<phi[37][57][35]<<"\n";
   std::cout.unsetf(std::ios_base::floatfield);
}

real frontPositionIn(real xyz1, real xyz2, real phi1, real phi2){
   real xyzfront;
   if (phi1 != phi2) {
      xyzfront=xyz1+(0.5-phi1)*(xyz2-xyz1)/(phi2-phi1);
   } else {
      xyzfront=0.5*(xyz1+xyz2);
      std::cout<<"#### phi1=phi2 ##############################\n";
   }
   return xyzfront;
}
#endif
