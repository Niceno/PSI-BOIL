#include "tif.h"

real Blending(const real tintnew, const real tifold, const real factor) {
  return  tintnew*factor + (tifold) * (1.0-factor);
}

/***************************************************************************//**
 *  Calculates interface temperature field
 *  Tn+1 = Tint*factor + Tn*(1-factor)
******************************************************************************/
void TIF::tint_field(const real factor, const bool iter) {
  if(variable_tif) {
    if(store_tif)
      for_vijk(tif,i,j,k) {
        if(!iter) {
          tifold[i][j][k] = tif[i][j][k];
        }
        tif   [i][j][k] = tr;
      }
    else
      for_vijk(tif,i,j,k) {
        tif   [i][j][k] = tr;
      }

    if(pres)
      Pressure_effect();

    Mass_src_effect();

    Extend_tint();

    if (store_tif) {
#if 1
      if(factor < 1.00) {
        for_vijk(tif,i,j,k)
          if(Interface(i,j,k) || Vicinity(i,j,k)) {
           real tintnew = tif[i][j][k];
            tif[i][j][k] = Blending(tintnew,tifold[i][j][k],factor);
          }
      }
#endif
    } else {
      for_vijk(tif,i,j,k) 
        tifold[i][j][k] = tif[i][j][k];
      boil::oout<<"EnthalpyFD::tint_field()  initialize tif"<<"\n";
      store_tif = true;
    } 
    tifold.exchange_all();
    tif   .exchange_all();
  }
}

void TIF::Extend_tint() { /* crude code */
  for_vijk(tif,i,j,k)
    if(!Interface(i,j,k) && Vicinity(i,j,k)) {
      int ctr(0);
      real tintf(0.0);

      if (Interface(i+1,j,k)) {
        ctr++;
        tintf += tif[i+1][j][k];
      } 
      if (Interface(i-1,j,k)) {
        ctr++;
        tintf += tif[i-1][j][k];
      } 
      if (Interface(i,j+1,k)) {
        ctr++;
        tintf += tif[i][j+1][k];
      } 
      if (Interface(i,j-1,k)) {
        ctr++;
        tintf += tif[i][j-1][k];
      } 
      if (Interface(i,j,k+1)) {
        ctr++;
        tintf += tif[i][j][k+1];
      } 
      if (Interface(i,j,k-1)) {
        ctr++;
        tintf += tif[i][j][k-1];
      }  

      if (ctr > 0) {
        tintf /= real(ctr);
        tif[i][j][k] = tintf;
      }
    }
}
