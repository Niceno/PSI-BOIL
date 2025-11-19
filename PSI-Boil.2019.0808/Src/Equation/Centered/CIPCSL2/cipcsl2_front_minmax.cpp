#include "cipcsl2.h"

real frontPosition(real xyz1, real xyz2, real phi1, real phi2);

/******************************************************************************/
void CIPCSL2::front_minmax() {
  front_minmax(Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa) );
}
/******************************************************************************/
void CIPCSL2::front_minmax(Range<real> xr,
                           Range<real> yr,
                           Range<real> zr ) {
/***************************************************************************//**
*  \brief Detect maximu and minimum of free surface position.
*           results: xminft,xmaxft,yminft,ymaxft,zminft,zmaxft
*******************************************************************************/

  xminft=boil::exa; xmaxft=-boil::exa;
  yminft=boil::exa; ymaxft=-boil::exa;
  zminft=boil::exa; zmaxft=-boil::exa;
  bool frontExist=false;

  real phi1, phi2, xyz1, xyz2, xfront, yfront, zfront;
  for (int i=phi.si()-1; i<=phi.ei()  ; i++) {
    if (phi.xc(i  )<xr.first()) continue;
    if (phi.xc(i+1)>xr.last() ) continue;
    for (int j=phi.sj()  ; j<=phi.ej()  ; j++) {
      if (phi.yc(j  )<yr.first()) continue;
      if (phi.yc(j  )>yr.last() ) continue;
      for (int k=phi.sk()  ; k<=phi.ek()  ; k++) {
        if (phi.zc(k  )<zr.first()) continue;
        if (phi.zc(k  )>zr.last() ) continue;
        /* i-direction */
        phi1=phi[i  ][j][k];
        phi2=phi[i+1][j][k];
#ifdef IB
        if ( dom->ibody().off(i,j,k) || dom->ibody().off(i+1,j,k)) continue;
#endif
        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           frontExist=true;
           xyz1=phi.xc(i);
           xyz2=phi.xc(i+1);
           xfront=frontPosition(xyz1,xyz2,phi1,phi2);
           yfront=phi.yc(j);
           zfront=phi.zc(k);
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) yminft=yfront;
           if(yfront>ymaxft) ymaxft=yfront;
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
        }
      }
    }
  }
  for (int i=phi.si()  ; i<=phi.ei()  ; i++) {
    if (phi.xc(i  )<xr.first()) continue;
    if (phi.xc(i  )>xr.last() ) continue;
    for (int j=phi.sj()-1; j<=phi.ej()  ; j++) {
      if (phi.yc(j  )<yr.first()) continue;
      if (phi.yc(j+1)>yr.last() ) continue;
      for (int k=phi.sk()  ; k<=phi.ek()  ; k++) {
        if (phi.zc(k  )<zr.first()) continue;
        if (phi.zc(k  )>zr.last() ) continue;
        /* j-direction */
        phi1=phi[i][j  ][k];
        phi2=phi[i][j+1][k];
#ifdef IB
        if ( dom->ibody().off(i,j,k) || dom->ibody().off(i,j+1,k)) continue;
#endif
        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           frontExist=true;
           xyz1=phi.yc(j);
           xyz2=phi.yc(j+1);
           xfront=phi.xc(i);
           yfront=frontPosition(xyz1,xyz2,phi1,phi2);
           zfront=phi.zc(k);
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) yminft=yfront;
           if(yfront>ymaxft) ymaxft=yfront;
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
        } 
      }
    }
  }
  for (int i=phi.si()  ; i<=phi.ei()  ; i++) {
    if (phi.xc(i  )<xr.first()) continue;
    if (phi.xc(i  )>xr.last() ) continue;
    for (int j=phi.sj()  ; j<=phi.ej()  ; j++) {
      if (phi.yc(j  )<yr.first()) continue;
      if (phi.yc(j  )>yr.last() ) continue;
      for (int k=phi.sk()-1; k<=phi.ek()  ; k++) {
        if (phi.zc(k  )<zr.first()) continue;
        if (phi.zc(k+1)>zr.last() ) continue;
        /* k-direction */
        phi1=phi[i][j][k  ];
        phi2=phi[i][j][k+1];
#ifdef IB
        if ( dom->ibody().off(i,j,k) || dom->ibody().off(i,j,k+1)) continue;
#endif
        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           frontExist=true;
           xyz1=phi.zc(k);
           xyz2=phi.zc(k+1);
           xfront=phi.xc(i);
           yfront=phi.yc(j);
           zfront=frontPosition(xyz1,xyz2,phi1,phi2);
           if(xfront<xminft) xminft=xfront;
           if(xfront>xmaxft) xmaxft=xfront;
           if(yfront<yminft) yminft=yfront;
           if(yfront>ymaxft) ymaxft=yfront;
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
   boil::oout << "cipcsl2_front_minmax:time,xmax= " << time->current_time() 
              << " " << xmaxft << " " << zmaxft << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);
}

real frontPosition(real xyz1, real xyz2, real phi1, real phi2){
   real xyzfront;
   if (phi1 != phi2) {
      xyzfront=xyz1+(0.5-phi1)*(xyz2-xyz1)/(phi2-phi1);
   } else {
      xyzfront=0.5*(xyz1+xyz2);
   } 
   return xyzfront;
}
