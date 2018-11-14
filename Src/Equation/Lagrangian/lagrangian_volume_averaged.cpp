#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::volume_averaged() {

  /*-----------------------------------------------+
  |  compute volume averaged velocity and density  |
  +-----------------------------------------------*/
  
  
  for_p(p) {
    for_m(m) { 
      real sum = 0.0; 
      real den = 0.0;
      real vol = 0.0;
      int ii=0; if(m == Comp::i()) ii++;
      int jj=0; if(m == Comp::j()) jj++;
      int kk=0; if(m == Comp::k()) kk++;
      for_pijk(p,i,j,k){
        const real c_a = interface_fraction(i,    j,    k,    p, continuous);
        const real c_b = interface_fraction(i-ii, j-jj, k-kk, p, continuous);
        const real c_av = 0.5*(c_a + c_b);
        const real d_frac = lagrangian_fraction(c_av);
        sum +=  (*u)[m][i][j][k]  * d_frac * (*u).dV(m,i,j,k); 
        den +=  flu->rho(m,i,j,k) * d_frac * (*u).dV(m,i,j,k);
        vol +=                      d_frac * (*u).dV(m,i,j,k);
      }
      boil::cart.sum_real(& sum);
      boil::cart.sum_real(& den);
      boil::cart.sum_real(& vol);
      particles[p].vol_rho(m, den/vol); 
      particles[p].vol_uvw(m, sum/vol); 
    }
  }

  OPR("lagrangian_volume_averaged.cpp");
  getchar();

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_volume_averaged.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
