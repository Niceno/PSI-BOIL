#include "nucleation.h"

real frontPos(real xyz1, real xyz2, real phi1, real phi2);

/******************************************************************************/
real Nucleation::zftmin( Range<real> xr
                       , Range<real> yr
                       , Range<real> zr ) {
/***************************************************************************//**
*  \brief Detect maximu and minimum of free surface position.
*           results: xminft,xmaxft,yminft,ymaxft,zminft,zmaxft
*******************************************************************************/

  real zminft=boil::exa;
  real zmaxft=-boil::exa;
  bool frontExist=false;
  bool bib=true;  // true: with IB,  false: without IB
  if ( clr->domain()->ibody().ncall()==0 )   bib=false;

  real phi1, phi2, xyz1, xyz2, xfront, yfront, zfront;
  for (int i=(*clr).si()  ; i<=(*clr).ei()  ; i++) {
    if ((*clr).xc(i  )<xr.first()) continue;
    if ((*clr).xc(i  )>xr.last() ) continue;
    for (int j=(*clr).sj()  ; j<=(*clr).ej()  ; j++) {
      if ((*clr).yc(j  )<yr.first()) continue;
      if ((*clr).yc(j  )>yr.last() ) continue;
      for (int k=(*clr).sk()-1; k<=(*clr).ek()  ; k++) {
        if ((*clr).zc(k+1)<zr.first()) continue;
        if ((*clr).zc(k  )>zr.last() ) continue;
        /* k-direction */
        phi1=(*clr)[i][j][k  ];
        phi2=(*clr)[i][j][k+1];

	/* exception: bubble touches to bottom */
        if ( clr->domain()->ibody().off(i,j,k) &&
             clr->domain()->ibody().on (i,j,k+1)) { // stride solid surface
          if(phi2<0.5) zminft=0.0;
        }
	if ( bib == false) {
          if ( k==1 && clr->zn(k)<zbtm+boil::pico ) {
            if(phi1<0.5) zminft=0.0;
          }
        }

	/* neglect solid region */
        if ( clr->domain()->ibody().off(i,j,k) ||
             clr->domain()->ibody().off(i,j,k+1)) continue;

        if ( (phi1-0.5)*(phi2-0.5) <=0.0) {
           frontExist=true;
           xyz1=(*clr).zc(k);
           xyz2=(*clr).zc(k+1);
           xfront=(*clr).xc(i);
           yfront=(*clr).yc(j);
           zfront=frontPos(xyz1,xyz2,phi1,phi2);
           if(zfront<zminft) zminft=zfront;
           if(zfront>zmaxft) zmaxft=zfront;
         }
       }
     }
   }
   boil::cart.min_real(&zminft);
   boil::cart.max_real(&zmaxft);
   return(zminft);
}

real frontPos(real xyz1, real xyz2, real phi1, real phi2){
   real xyzfront;
   if (phi1 != phi2) {
      xyzfront=xyz1+(0.5-phi1)*(xyz2-xyz1)/(phi2-phi1);
   } else {
      xyzfront=0.5*(xyz1+xyz2);
   } 
   return xyzfront;
}
/*-----------------------------------------------------------------------------+
 '$Id: nucleation_zftmin.cpp,v 1.1 2014/08/06 08:19:36 sato Exp $'/
+-----------------------------------------------------------------------------*/
