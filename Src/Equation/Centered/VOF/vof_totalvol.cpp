#include "vof.h"

/******************************************************************************/
real VOF::totalvol(real * vapss) {
  return totalvol(Range<real>(-boil::exa, boil::exa),
                  Range<real>(-boil::exa, boil::exa),
                  Range<real>(-boil::exa, boil::exa),
                  vapss);
}
/******************************************************************************/
real VOF::totalvol(Range<real> xr, Range<real> yr, Range<real> zr,
                   real * vapss) {
   real phisum = 0.0;
   real vapsum = 0.0;

   for_ijk(i,j,k){
     if(dom->ibody().off(i,j,k)) continue;
     if (phi.xn(i  )<xr.first()) continue;
     if (phi.xn(i+1)>xr.last() ) continue;
     if (phi.yn(j  )<yr.first()) continue;
     if (phi.yn(j+1)>yr.last() ) continue;
     if (phi.zn(k  )<zr.first()) continue;
     if (phi.zn(k+1)>zr.last() ) continue;
     phisum += phi[i][j][k] * dV(i,j,k);
     vapsum += (1.0-phi[i][j][k]) * dV(i,j,k);
   }

   boil::cart.sum_real(&phisum);
   boil::cart.sum_real(&vapsum);


   std::cout.setf(std::ios_base::scientific);
   boil::oout << "totalvol:time,volume,phisum= " 
              << time->current_time()
              <<" "<< phisum << " "<<vapsum<< boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);

   if(vapss)
     *vapss = vapsum;

   return phisum;
}
