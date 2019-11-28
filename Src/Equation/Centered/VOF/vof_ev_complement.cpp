#include "vof.h"

/******************************************************************************/
void VOF::ev_complement(const ScalarInt & pflag, const Scalar & scp,
                        const Vector & umixed, const Vector & uliq, 
                        Vector & ugas, const Vector * bndclr) {
/***************************************************************************//**
*  \brief Calculate the velocity correction for the complementary phase.
*******************************************************************************/

  for_m(m)
    for_avmijk(umixed,m,i,j,k)
      ugas[m][i][j][k] = umixed[m][i][j][k];

  Comp m;

  m = Comp::u();
  if(ifull)
    for_vmijk(ugas,m,i,j,k) {
      if(uliq[m][i][j][k]!=umixed[m][i][j][k]) {
        real bdval;
        if(bndclr) {
          bdval = (*bndclr)[m][i][j][k];
        } else {
          bdval = std::min(1.0,std::max(0.0,0.5*(scp[i-1][j][k]+scp[i][j][k])));
        }
        if((1.-bdval)>boil::pico) {
          ugas[m][i][j][k] = (umixed[m][i][j][k]-bdval*uliq[m][i][j][k])
                             /(1.-bdval);
        }
      }
    }

  m = Comp::v();
  if(jfull)
    for_vmijk(ugas,m,i,j,k) {
      if(uliq[m][i][j][k]!=umixed[m][i][j][k]) {
        real bdval;
        if(bndclr) {
          bdval = (*bndclr)[m][i][j][k];
        } else {
          bdval = std::min(1.0,std::max(0.0,0.5*(scp[i][j-1][k]+scp[i][j][k])));
        }
        if((1.-bdval)>boil::pico) {
          ugas[m][i][j][k] = (umixed[m][i][j][k]-bdval*uliq[m][i][j][k])
                             /(1.-bdval);
        }              
      }
    }

  m = Comp::w();
  if(kfull)
    for_vmijk(ugas,m,i,j,k) {
      if(uliq[m][i][j][k]!=umixed[m][i][j][k]) {
        real bdval;
        if(bndclr) {
          bdval = (*bndclr)[m][i][j][k];
        } else {
          bdval = std::min(1.0,std::max(0.0,0.5*(scp[i][j][k-1]+scp[i][j][k])));
        }
        if((1.-bdval)>boil::pico) {
          ugas[m][i][j][k] = (umixed[m][i][j][k]-bdval*uliq[m][i][j][k])
                             /(1.-bdval);
        }
      }
    }

  ugas.exchange_all();

  return;
}
