#include "topology.h"

/******************************************************************************/
void Topology::front_minmax(real * store_arr) {
  front_minmax(Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               store_arr);
}

/******************************************************************************/
void Topology::front_minmax(Range<real> xr,
                            Range<real> yr,
                            Range<real> zr,
                            real * store_arr) {
/***************************************************************************//**
*  \brief Detect maximum and minimum of free surface position.
*         results: xminft_tmp,xmaxft_tmp,yminft_tmp,ymaxft_tmp,zminft_tmp,zmaxft_tmp
*******************************************************************************/

   real xminft_tmp=boil::exa; real xmaxft_tmp=-boil::exa;
   real yminft_tmp=boil::exa; real ymaxft_tmp=-boil::exa;
   real zminft_tmp=boil::exa; real zmaxft_tmp=-boil::exa;
   bool frontExist=false;

   real xfront, yfront, zfront;
   real phim, phip;

   /* i-direction */
   for(int i=clr->si(); i<=clr->ei()+1; i++) {
     if ((*clr).xc(i  )<xr.first()) continue;
     if ((*clr).xc(i+1)>xr.last() ) continue;
     for(int j=clr->sj(); j<=clr->ej()  ; j++) {
       if ((*clr).yc(j  )<yr.first()) continue;
       if ((*clr).yc(j  )>yr.last() ) continue;
       for(int k=clr->sk(); k<=clr->ek()  ; k++) {
         if ((*clr).zc(k  )<zr.first()) continue;
         if ((*clr).zc(k  )>zr.last() ) continue;
         if ( clr->domain()->ibody().off(i,j,k) ||
              clr->domain()->ibody().off(i+1,j,k)) continue;
         phim=(*clr)[i-1][j][k];
         phip=(*clr)[i  ][j][k];
         if((phim-clrsurf)*(phip-clrsurf)<=0.0) {
           frontExist=true;
           xfront=frontPosition(i,j,k,Comp::i());
           yfront=(*clr).yc(j);
           zfront=(*clr).zc(k);
           if(xfront<xminft_tmp) xminft_tmp=xfront;
           if(xfront>xmaxft_tmp) xmaxft_tmp=xfront;
           if(yfront<yminft_tmp) yminft_tmp=yfront;
           if(yfront>ymaxft_tmp) ymaxft_tmp=yfront;
           if(zfront<zminft_tmp) zminft_tmp=zfront;
           if(zfront>zmaxft_tmp) zmaxft_tmp=zfront;
         }
       }  /* k */
     }  /* j */
   } /* i */

   /* j-direction */
   for(int i=clr->si(); i<=clr->ei()  ; i++) {
     if ((*clr).xc(i  )<xr.first()) continue;
     if ((*clr).xc(i  )>xr.last() ) continue;
     for(int j=clr->sj(); j<=clr->ej()+1; j++) {
       if ((*clr).yc(j  )<yr.first()) continue;
       if ((*clr).yc(j+1)>yr.last() ) continue;
       for(int k=clr->sk(); k<=clr->ek()  ; k++) {
         if ((*clr).zc(k  )<zr.first()) continue;
         if ((*clr).zc(k  )>zr.last() ) continue;
         if ( clr->domain()->ibody().off(i,j,k) ||
              clr->domain()->ibody().off(i,j+1,k)) continue;
         phim=(*clr)[i][j-1][k];
         phip=(*clr)[i][j  ][k];
         if((phim-clrsurf)*(phip-clrsurf)<=0.0) {
           frontExist=true;
           xfront=(*clr).xc(i);
           yfront=frontPosition(i,j,k,Comp::j());
           zfront=(*clr).zc(k);
           if(xfront<xminft_tmp) xminft_tmp=xfront;
           if(xfront>xmaxft_tmp) xmaxft_tmp=xfront;
           if(yfront<yminft_tmp) yminft_tmp=yfront;
           if(yfront>ymaxft_tmp) ymaxft_tmp=yfront;
           if(zfront<zminft_tmp) zminft_tmp=zfront;
           if(zfront>zmaxft_tmp) zmaxft_tmp=zfront;
         } 
       } /* k */
     } /* j */
   } /* i */

   /* k-direction */
   for(int i=clr->si(); i<=clr->ei()  ; i++) {
     if ((*clr).xc(i  )<xr.first()) continue;
     if ((*clr).xc(i  )>xr.last() ) continue;
     for(int j=clr->sj(); j<=clr->ej()  ; j++) {
       if ((*clr).yc(j  )<yr.first()) continue;
       if ((*clr).yc(j  )>yr.last() ) continue;
       for(int k=clr->sk(); k<=clr->ek()+1; k++) {
         if ((*clr).zc(k  )<zr.first()) continue;
         if ((*clr).zc(k+1)>zr.last() ) continue;
         if ( clr->domain()->ibody().off(i,j,k) ||
              clr->domain()->ibody().off(i,j,k+1)) continue;
         phim=(*clr)[i][j][k-1];
         phip=(*clr)[i][j][k  ];
         if((phim-clrsurf)*(phip-clrsurf)<=0.0) {
           frontExist=true;
           xfront=(*clr).xc(i);
           yfront=(*clr).yc(j);
           zfront=frontPosition(i,j,k,Comp::k());
           if(xfront<xminft_tmp) xminft_tmp=xfront;
           if(xfront>xmaxft_tmp) xmaxft_tmp=xfront;
           if(yfront<yminft_tmp) yminft_tmp=yfront;
           if(yfront>ymaxft_tmp) ymaxft_tmp=yfront;
           if(zfront<zminft_tmp) zminft_tmp=zfront;
           if(zfront>zmaxft_tmp) zmaxft_tmp=zfront;
         }
       } /* k */
     } /* j */
   } /* i */

   boil::cart.min_real(&xminft_tmp);
   boil::cart.max_real(&xmaxft_tmp);
   boil::cart.min_real(&yminft_tmp);
   boil::cart.max_real(&ymaxft_tmp);
   boil::cart.min_real(&zminft_tmp);
   boil::cart.max_real(&zmaxft_tmp);

   if(store_arr) {
     store_arr[0] = xminft_tmp;
     store_arr[1] = xmaxft_tmp;
     store_arr[2] = yminft_tmp;
     store_arr[3] = ymaxft_tmp;
     store_arr[4] = zminft_tmp;
     store_arr[5] = zmaxft_tmp;
   } else {
     xminft = xminft_tmp;
     xmaxft = xmaxft_tmp;
     yminft = yminft_tmp;
     ymaxft = ymaxft_tmp;
     zminft = zminft_tmp;
     zmaxft = zmaxft_tmp;
   }

   return;
}

real Topology::frontPosition(const int i, const int j, const int k,
                             const Comp m) {
   real xyzfront;
   if (m == Comp::i()) { 
     xyzfront = (*fs)[m][i][j][k];
     if(!boil::realistic(xyzfront))
       xyzfront = (*clr).xn(i);
   } else if (m == Comp::j()) {
     xyzfront = (*fs)[m][i][j][k];
     if(!boil::realistic(xyzfront))
       xyzfront = (*clr).yn(j);
   } else {
     xyzfront = (*fs)[m][i][j][k];
     if(!boil::realistic(xyzfront))
       xyzfront = (*clr).zn(k);
   }

   return xyzfront;
}
