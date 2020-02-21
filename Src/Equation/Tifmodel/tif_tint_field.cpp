#include "tif.h"

inline real TIF::underrelaxation(const real tintnew, const real tifold) {
  return  tintnew*factor + (tifold) * (1.0-factor);
}

/***************************************************************************//**
 *  Calculates interface temperature field
 *  Tn+1 = Tint*factor + Tn*(1-factor)
*******************************************************************************/
void TIF::tint_field(const bool newstep) {
  if(variable_tif) {
    if(store_tif) {
      if(newstep) {
        tifold = tif;
      }
    }

    /* individual tifmodel */
    model();

    if(weaklim) {
      for_vijk(tif,i,j,k) {
        if(Interface(i,j,k)) {
          tif[i][j][k] = std::min(tmax,
                                  std::max(tmin,tif[i][j][k])
                                  );
        }
      }
    }

    tif.bnd_update();
    tif.exchange_all();

    /* under-relaxation */
    if(store_tif) {
      if(factor < 1.00) {
        for_vijk(tif,i,j,k) {
          if(Interface(i,j,k)) {
            real tintnew = tif[i][j][k];
            tif[i][j][k] = underrelaxation(tintnew,tifold[i][j][k]);
          }
        }
      }
    } 
   
    tif.bnd_update();
    tif.exchange_all();

    /* extend the field */
    extend_tint();
  
    if(!store_tif) {
      tifold = tif;
      boil::oout<<"TIFmodel::tint_field()  initialize tif"<<"\n";
      store_tif = true;
    } 

  } /* variable */

  return;
}
