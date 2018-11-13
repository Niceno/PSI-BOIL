#include "enthalpytif.h"

real Blending(const real tintnew, const real tifold, const real factor) {
  return  tintnew*factor + (tifold) * (1.0-factor);
}

/***************************************************************************//**
 *  Calculates interface temperature field
 *  Tn+1 = Tint*factor + Tn*(1-factor)
******************************************************************************/
void EnthalpyTIF::tint_field(const Scalar & heaviside, const real factor, const bool iter) {
  bool trigger = false;
  if(pres||mflx)
    trigger = true;

  if(trigger) {
    if(store_tif)
      for_ijk(i,j,k) {
        if(!iter) {
          tifold[i][j][k] = tif[i][j][k];
        }
        tif   [i][j][k] = tsat;
      }
    else
      for_ijk(i,j,k) {
        tif   [i][j][k] = tsat;
      }

    if(pres)
      Pressure_effect(heaviside);
    if(mflx)
      Mass_src_effect(heaviside);

    Extend_tint(heaviside);

    if (store_tif) {
#if 1
      if(factor < 1.00) {
        for_ijk(i,j,k)
          if(Interface(i,j,k,heaviside) || Vicinity(i,j,k,heaviside)) {
           real tintnew = tif[i][j][k];
            tif[i][j][k] = Blending(tintnew,tifold[i][j][k],factor);
          }
      }
#endif
    } else {
      for_ijk(i,j,k) 
        tifold[i][j][k] = tif[i][j][k];
      boil::oout<<"EnthalpyFD::tint_field()  initialize tif"<<"\n";
      store_tif = true;
    } 
    tifold.exchange_all();
    tif   .exchange_all();
  }
}

void EnthalpyTIF::Extend_tint(const Scalar & heaviside) {
  for_ijk(i,j,k)
    if(!Interface(i,j,k,heaviside) && Vicinity(i,j,k,heaviside)) {
      int ctr(0);
      real tintf(0.0);

      if (Interface(heaviside[i+1][j][k])) {
        ctr++;
        tintf += tif[i+1][j][k];
      } 
      if (Interface(heaviside[i-1][j][k])) {
        ctr++;
        tintf += tif[i-1][j][k];
      } 
      if (Interface(heaviside[i][j+1][k])) {
        ctr++;
        tintf += tif[i][j+1][k];
      } 
      if (Interface(heaviside[i][j-1][k])) {
        ctr++;
        tintf += tif[i][j-1][k];
      } 
      if (Interface(heaviside[i][j][k+1])) {
        ctr++;
        tintf += tif[i][j][k+1];
      } 
      if (Interface(heaviside[i][j][k-1])) {
        ctr++;
        tintf += tif[i][j][k-1];
      }  

      if (ctr > 0) {
        tintf /= real(ctr);
        tif[i][j][k] = tintf;
      }
    }
}
