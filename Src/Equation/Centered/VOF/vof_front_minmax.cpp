#include "vof.h"

real frontPosition(real xyz1, real xyz2, real phi1, real phi2);

/******************************************************************************/
void VOF::front_minmax() {
/***************************************************************************//**
*  \brief Detect maximu and minimum of free surface position.
*           results: xminft,xmaxft,yminft,ymaxft,zminft,zmaxft
*******************************************************************************/

   xminft=1.0E300; xmaxft=-1.0E300;
   yminft=1.0E300; ymaxft=-1.0E300;
   zminft=1.0E300; zmaxft=-1.0E300;
   bool frontExist=false;

   real phi1, phi2, xyz1, xyz2, xfront, yfront, zfront;
   for (int i=phi.si()-1; i<=phi.ei()-1; i++) {
      for (int j=phi.sj()-1; j<=phi.ej()-1; j++) {
         for (int k=phi.sk()-1; k<=phi.ek()-1; k++) {
            /* i-direction */
               phi1=phi[i  ][j][k];
               phi2=phi[i+1][j][k];
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
            /* j-direction */
               phi1=phi[i][j  ][k];
               phi2=phi[i][j+1][k];
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
            /* k-direction */
               phi1=phi[i][j][k  ];
               phi2=phi[i][j][k+1];
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
   boil::oout << "front_minmax:time,xmax= " << time->current_time() 
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
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_front_minmax.cpp,v 1.3 2009/11/12 12:15:48 sato Exp $'/
+-----------------------------------------------------------------------------*/
