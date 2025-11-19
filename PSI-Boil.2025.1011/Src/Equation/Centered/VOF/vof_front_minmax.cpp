#include "vof.h"

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
