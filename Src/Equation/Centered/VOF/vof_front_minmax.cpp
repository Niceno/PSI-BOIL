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

   real xfront, yfront, zfront;
   real phim, phip;

   /* i-direction */
   for(int i=si(); i<=ei()+1; i++)
   for(int j=sj(); j<=ej()  ; j++)
   for(int k=sk(); k<=ek()  ; k++) {

     phim=phi[i-1][j][k];
     phip=phi[i  ][j][k];
     if((phim-phisurf)*(phip-phisurf)<=0.0) {
       frontExist=true;
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
   }

   /* j-direction */
   for(int i=si(); i<=ei()  ; i++)
   for(int j=sj(); j<=ej()+1; j++)
   for(int k=sk(); k<=ek()  ; k++) {
  
     phim=phi[i][j-1][k];
     phip=phi[i][j  ][k];
     if((phim-phisurf)*(phip-phisurf)<=0.0) {
       frontExist=true;
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
   }

   /* k-direction */
   for(int i=si(); i<=ei()  ; i++)
   for(int j=sj(); j<=ej()  ; j++)
   for(int k=sk(); k<=ek()+1; k++) {
  
     phim=phi[i][j][k-1];
     phip=phi[i][j][k  ];
     if((phim-phisurf)*(phip-phisurf)<=0.0) {
       frontExist=true;
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

real VOF::frontPosition(const int i, const int j, const int k,
                        const Comp m){
   real xyzfront;
   if (m == Comp::i()) { 
     xyzfront = fs[m][i][j][k];
     if (xyzfront>boil::zetta)
       xyzfront = phi.xn(i);
   } else if (m == Comp::j()) {
     xyzfront = fs[m][i][j][k];
     if (xyzfront>boil::zetta)
       xyzfront = phi.yn(j);
   } else {
     xyzfront = fs[m][i][j][k];
     if (xyzfront>boil::zetta)
       xyzfront = phi.zn(k);
   }

   return xyzfront;
}
