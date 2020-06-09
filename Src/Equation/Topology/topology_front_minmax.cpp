#include "topology.h"

/******************************************************************************/
void Topology::front_minmax() {
  front_minmax(Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa) );
}

/******************************************************************************/
void Topology::front_minmax(Range<real> xr,
                            Range<real> yr,
                            Range<real> zr ) {
/***************************************************************************//**
*  \brief Detect maximum and minimum of free surface position.
*         results: xminft,xmaxft,yminft,ymaxft,zminft,zmaxft
*******************************************************************************/

   xminft=boil::exa; xmaxft=-boil::exa;
   yminft=boil::exa; ymaxft=-boil::exa;
   zminft=boil::exa; zmaxft=-boil::exa;
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
