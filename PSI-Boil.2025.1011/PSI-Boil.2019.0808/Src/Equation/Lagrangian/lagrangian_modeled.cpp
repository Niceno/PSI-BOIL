#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::modeled(Vector * force) {
/*-----------------------------------------------+
|   compute the correcting forces which will be  |
|    added as a source term in the momentum eq   |
+-----------------------------------------------*/

  for_p(p) { 
    real duvw = 0.0;

    for_m(m) {

      duvw = particles[p].uvw(m) - particles[p].vol_uvw(m);
      int si = particles[p].si(); int ei = particles[p].ei();  
      int sj = particles[p].sj(); int ej = particles[p].ej();  
      int sk = particles[p].sk(); int ek = particles[p].ek();  
      int ii=0; int jj=0; int kk=0; 
      if(m == Comp::i()) {ei++; ii++;}
      if(m == Comp::j()) {ej++; jj++;}
      if(m == Comp::k()) {ek++; kk++;}

      for(int i=si; i<=ei; i++)
        for(int j=sj; j<=ej; j++)
          for(int k=sk; k<=ek; k++) {

            const real c_a = interface_fraction(i,j,k,p,continuous);
            const real c_b = interface_fraction(i-ii,j-jj,k-kk,p,continuous);
            const real c_av = 0.5*(c_a + c_b);
            const real d_frac = lagrangian_fraction(c_av);
            if(d_frac > 0.5) { 
               (*force)[m][i][j][k] += 
                 particles[p].vol_rho(m) * duvw / time->dt()
                                         * force->dV(m,i,j,k)
                                         * d_frac ;
            }
      } // for i,j,k
    } // for_m
  } // for_p

  OPR("lagrangian_modeled.cpp");
  getchar();

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_modeled.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
