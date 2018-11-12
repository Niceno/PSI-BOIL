#include "cipcsl2.h"
#include <iomanip>

/******************************************************************************/
void CIPCSL2::totalvol() {
/***************************************************************************//**
*  \brief Calculate volume by two methods.
*         Method1: Sum(clr*dV)
*         Method2: Sum(dV)  in the region clr > phisurf
*******************************************************************************/

   clrsum1=0.0;
   clrsum2=0.0;

   if(dom->ibody().nccells()==0){
     for_ijk(i,j,k){
       clrsum1 += clr[i][j][k] * dV(i,j,k);
       clrsum2 += (1.0-clr[i][j][k]) * dV(i,j,k);
     }
   } else {
     for_ijk(i,j,k){
       clrsum1 += clr[i][j][k]       * dV(i,j,k) * dom->ibody().fV(i,j,k);
       clrsum2 += (1.0-clr[i][j][k]) * dV(i,j,k) * dom->ibody().fV(i,j,k);
     }
   }

   boil::cart.sum_real(&clrsum1);
   boil::cart.sum_real(&clrsum2);

#if 0
   std::cout.setf(std::ios_base::scientific);
   std::cout<< std::setprecision(16);
   boil::oout << "cipcsl2_totalvol:time,clrsum1,clrsum2= " 
              << time->current_time()
              <<" "<< clrsum1 <<" "<< clrsum2 << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);
   std::cout<< std::setprecision(6);
#endif

   return;
}

/******************************************************************************/
void CIPCSL2::totalvol( Range<real> xr, Range<real> yr, Range<real> zr) {
/***************************************************************************//**
*  \brief Calculate volume by two methods.
*         Method1: Sum(clr*dV)
*         Method2: Sum(dV)  in the region clr > phisurf
*******************************************************************************/

   real clrsum1r=0.0;
   real clrsum2r=0.0;

   if(dom->ibody().nccells()==0){
     for (int i=phi.si(); i<=phi.ei()  ; i++) {
       if (phi.xc(i)<xr.first()) continue;
       if (phi.xc(i)>xr.last() ) continue;
       for (int j=phi.sj()  ; j<=phi.ej()  ; j++) {
         if (phi.yc(j)<yr.first()) continue;
         if (phi.yc(j)>yr.last() ) continue;
         for (int k=phi.sk()  ; k<=phi.ek()  ; k++) {
           if (phi.zc(k)<zr.first()) continue;
           if (phi.zc(k)>zr.last() ) continue;
           clrsum1r += clr[i][j][k] * dV(i,j,k);
           clrsum2r += (1.0-clr[i][j][k]) * dV(i,j,k);
         }
       }
     }
   } else {
     for (int i=phi.si(); i<=phi.ei()  ; i++) {
       if (phi.xc(i)<xr.first()) continue;
       if (phi.xc(i)>xr.last() ) continue;
       for (int j=phi.sj()  ; j<=phi.ej()  ; j++) {
         if (phi.yc(j)<yr.first()) continue;
         if (phi.yc(j)>yr.last() ) continue;
         for (int k=phi.sk()  ; k<=phi.ek()  ; k++) {
           if (phi.zc(k)<zr.first()) continue;
           if (phi.zc(k)>zr.last() ) continue;
           clrsum1r += clr[i][j][k]       * dV(i,j,k) * dom->ibody().fV(i,j,k);
           clrsum2r += (1.0-clr[i][j][k]) * dV(i,j,k) * dom->ibody().fV(i,j,k);
         }
       }
     }
   }

   boil::cart.sum_real(&clrsum1r);
   boil::cart.sum_real(&clrsum2r);

   boil::oout << "cipcsl2_totalvolRange:time,clrsum1,clrsum2= " 
              << time->current_time()
              <<" "<< clrsum1r <<" "<< clrsum2r <<"\n";

   return;
}




/******************************************************************************/
real CIPCSL2::totalvol(const Scalar & sca ) {
/***************************************************************************//**
*  \brief Calculate volume.
*         input:sca
*******************************************************************************/

   real scasum = 0.0;

   for_ijk(i,j,k)
     scasum += sca[i][j][k] * dV(i,j,k);

   boil::cart.sum_real(&scasum);

   return(scasum);
}
