#include "momentum.h"
#include "../../../Plot/plot.h"

/******************************************************************************/
void Momentum::convection(Vector * conv) {
/*------------+
|       3  n  |
|  f += - H   |
|       2     |
+------------*/

  Comp m;

  u.exchange();

  for_m(m)
    for(int i=0; i<u.ni(m); i++)
      for(int j=0; j<u.nj(m); j++)
        for(int k=0; k<u.nk(m); k++) {
          (*conv)[m][i][j][k] = 0.0;
            fnew [m][i][j][k] = 0.0;
        }

  /////////
  //     //
  //  u  //
  //     //
  /////////
  m = Comp::u();

  /////////////
  //// u u ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* u */
    real a_w = dSx(m,Sign::neg(),i,j,k); 
    real a_e = dSx(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_w *= dom->ibody().fSw(m,i,j,k);
      a_e *= dom->ibody().fSe(m,i,j,k);
    } 

    real umf = 0.5 * (u[Comp::u()][i-1][j][k] + u[Comp::u()][i][j][k]);      // u @ imin
    real upf = 0.5 * (u[Comp::u()][i+1][j][k] + u[Comp::u()][i][j][k]);      // u @ imax
    
//    um = lim.limit(-umf, u[Comp::u()][i+1][j][k],u[Comp::u()][i][j][k],u[Comp::u()][i-1][j][k]); 
//    up = lim.limit(+upf, u[Comp::u()][i-1][j][k],u[Comp::u()][i][j][k],u[Comp::u()][i+1][j][k]); 
    real um = 0.5 * (u[Comp::u()][i][j][k]+u[Comp::u()][i-1][j][k]);
    real up = 0.5 * (u[Comp::u()][i][j][k]+u[Comp::u()][i+1][j][k]);

    (*conv)[Comp::u()][i]  [j][k] += a_w * umf*um - a_e * upf*up;
    } 

  /////////////
  //// v u ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* v */
    real a_s = dSy(m,Sign::neg(),i,j,k); 
    real a_n = dSy(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_s *= dom->ibody().fSs(m,i,j,k);
      a_n *= dom->ibody().fSn(m,i,j,k);
    } 

    real vmf = 0.5 * (u[Comp::v()][i-1][j]  [k] + u[Comp::v()][i][j]  [k]);  // v @ jmin
    real vpf = 0.5 * (u[Comp::v()][i-1][j+1][k] + u[Comp::v()][i][j+1][k]);  // v @ jmax

//    um = lim.limit(-vmf, u[Comp::u()][i][j+1][k],u[Comp::u()][i][j][k],u[Comp::u()][i][j-1][k]); 
//    up = lim.limit(+vpf, u[Comp::u()][i][j-1][k],u[Comp::u()][i][j][k],u[Comp::u()][i][j+1][k]);
    real um = 0.5 * (u[Comp::u()][i][j][k]+u[Comp::u()][i][j-1][k]);
    real up = 0.5 * (u[Comp::u()][i][j][k]+u[Comp::u()][i][j+1][k]);
    
    (*conv)[Comp::u()][i][j]  [k] += a_s * vmf*um - a_n * vpf*up;
    } 

  /////////////
  //// w u ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* w */
    real a_b = dSz(m,Sign::neg(),i,j,k); 
    real a_t = dSz(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_b *= dom->ibody().fSb(m,i,j,k);
      a_t *= dom->ibody().fSt(m,i,j,k);
    } 

    real wmf = 0.5 * (u[Comp::w()][i-1][j][k]   + u[Comp::w()][i][j][k]);    // w @ kmin
    real wpf = 0.5 * (u[Comp::w()][i-1][j][k+1] + u[Comp::w()][i][j][k+1]);  // w @ kmax

//    um = lim.limit(-wmf, u[Comp::u()][i][j][k+1],u[Comp::u()][i][j][k],u[Comp::u()][i][j][k-1]);
//    up = lim.limit(+wpf, u[Comp::u()][i][j][k-1],u[Comp::u()][i][j][k],u[Comp::u()][i][j][k+1]);
    real um = 0.5 * (u[Comp::u()][i][j][k]+u[Comp::u()][i][j][k-1]);
    real up = 0.5 * (u[Comp::u()][i][j][k]+u[Comp::u()][i][j][k+1]);

    (*conv)[Comp::u()][i][j][k]   += a_b * wmf*um - a_t * wpf*up;
    }
  
  /////////
  //     //
  //  v  //
  //     //
  /////////
  m = Comp::v();

  /////////////
  //// v v ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* v */
    real a_s = dSy(m,Sign::neg(),i,j,k); 
    real a_n = dSy(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_s *= dom->ibody().fSs(m,i,j,k);
      a_n *= dom->ibody().fSn(m,i,j,k);
    } 

    real vmf = 0.5 * (u[Comp::v()][i][j-1][k] + u[Comp::v()][i][j][k]);
    real vpf = 0.5 * (u[Comp::v()][i][j+1][k] + u[Comp::v()][i][j][k]);
    
//    vm = lim.limit(-vmf, u[Comp::v()][i][j+1][k],u[Comp::v()][i][j][k],u[Comp::v()][i][j-1][k]); 
//    vp = lim.limit(+vpf, u[Comp::v()][i][j-1][k],u[Comp::v()][i][j][k],u[Comp::v()][i][j+1][k]);
    real vm = 0.5 * (u[Comp::v()][i][j][k]+u[Comp::v()][i][j-1][k]);
    real vp = 0.5 * (u[Comp::v()][i][j][k]+u[Comp::v()][i][j+1][k]);

    (*conv)[Comp::v()][i][j]  [k] += a_s * vmf*vm - a_n * vpf*vp;
    }

  /////////////
  //// u v ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* u */
    real a_w = dSx(m,Sign::neg(),i,j,k); 
    real a_e = dSx(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_w *= dom->ibody().fSw(m,i,j,k);
      a_e *= dom->ibody().fSe(m,i,j,k);
    } 

    real umf = 0.5 * (u[Comp::u()][i]  [j][k] + u[Comp::u()][i]  [j-1][k]);  // u @ imin
    real upf = 0.5 * (u[Comp::u()][i+1][j][k] + u[Comp::u()][i+1][j-1][k]);  // u @ imax

//    vm = lim.limit(-umf, u[Comp::v()][i+1][j][k],u[Comp::v()][i][j][k],u[Comp::v()][i-1][j][k]);
//    vp = lim.limit(+upf, u[Comp::v()][i-1][j][k],u[Comp::v()][i][j][k],u[Comp::v()][i+1][j][k]);
    real vm = 0.5 * (u[Comp::v()][i][j][k]+u[Comp::v()][i-1][j][k]);
    real vp = 0.5 * (u[Comp::v()][i][j][k]+u[Comp::v()][i+1][j][k]);

    (*conv)[Comp::v()][i]  [j][k] += a_w * umf*vm - a_e * upf*vp;
    }

  /////////////
  //// w v ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* w */
    real a_b = dSz(m,Sign::neg(),i,j,k); 
    real a_t = dSz(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_b *= dom->ibody().fSb(m,i,j,k);
      a_t *= dom->ibody().fSt(m,i,j,k);
    } 

    real wmf = 0.5 * (u[Comp::w()][i][j][k]   + u[Comp::w()][i][j-1][k]);    // w @ kmin
    real wpf = 0.5 * (u[Comp::w()][i][j][k+1] + u[Comp::w()][i][j-1][k+1]);  // w @ kmax

//    vm = lim.limit(-wmf, u[Comp::v()][i][j][k+1],u[Comp::v()][i][j][k],u[Comp::v()][i][j][k-1]); 
//    vp = lim.limit(+wpf, u[Comp::v()][i][j][k-1],u[Comp::v()][i][j][k],u[Comp::v()][i][j][k+1]);
    real vm = 0.5 * (u[Comp::v()][i][j][k]+u[Comp::v()][i][j][k-1]);
    real vp = 0.5 * (u[Comp::v()][i][j][k]+u[Comp::v()][i][j][k+1]);

    (*conv)[Comp::v()][i][j][k]   += a_b * wmf*vm - a_t * wpf*vp;
    }

  /////////
  //     //
  //  w  //
  //     //
  /////////
  m = Comp::w();

  /////////////
  //// w w ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* w */
    real a_b = dSz(m,Sign::neg(),i,j,k); 
    real a_t = dSz(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_b *= dom->ibody().fSb(m,i,j,k);
      a_t *= dom->ibody().fSt(m,i,j,k);
    } 

    real wmf = 0.5 * (u[Comp::w()][i][j][k-1] + u[Comp::w()][i][j][k]);      // w @ kmin
    real wpf = 0.5 * (u[Comp::w()][i][j][k+1] + u[Comp::w()][i][j][k]);      // w @ kmax
    
//    wm = lim.limit(-wmf, u[Comp::w()][i][j][k+1],u[Comp::w()][i][j][k],u[Comp::w()][i][j][k-1]); 
//    wp = lim.limit(+wpf, u[Comp::w()][i][j][k-1],u[Comp::w()][i][j][k],u[Comp::w()][i][j][k+1]);
    real wm = 0.5 * (u[Comp::w()][i][j][k]+u[Comp::w()][i][j][k-1]);
    real wp = 0.5 * (u[Comp::w()][i][j][k]+u[Comp::w()][i][j][k+1]);

    (*conv)[Comp::w()][i][j][k]   += a_b * wmf*wm - a_t * wpf*wp;
    }

  /////////////
  //// u w ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* u */
    real a_w = dSx(m,Sign::neg(),i,j,k); 
    real a_e = dSx(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_w *= dom->ibody().fSw(m,i,j,k);
      a_e *= dom->ibody().fSe(m,i,j,k);
    } 

    real umf = 0.5 * (u[Comp::u()][i]  [j][k-1] + u[Comp::u()][i]  [j][k]);  // u @ imin
    real upf = 0.5 * (u[Comp::u()][i+1][j][k-1] + u[Comp::u()][i+1][j][k]);  // u @ imax

//    wm = lim.limit(-umf, u[Comp::w()][i+1][j][k],u[Comp::w()][i][j][k],u[Comp::w()][i-1][j][k]);
//    wp = lim.limit(+upf, u[Comp::w()][i-1][j][k],u[Comp::w()][i][j][k],u[Comp::w()][i+1][j][k]);
    real wm = 0.5 * (u[Comp::w()][i][j][k]+u[Comp::w()][i-1][j][k]);
    real wp = 0.5 * (u[Comp::w()][i][j][k]+u[Comp::w()][i+1][j][k]);

    (*conv)[Comp::w()][i]  [j][k] += a_w * umf*wm - a_e * upf*wp;
    }

  /////////////
  //// v w ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* v */
    real a_s = dSy(m,Sign::neg(),i,j,k); 
    real a_n = dSy(m,Sign::pos(),i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_s *= dom->ibody().fSs(m,i,j,k);
      a_n *= dom->ibody().fSn(m,i,j,k);
    } 

    real vmf = 0.5 * (u[Comp::v()][i][j]  [k-1] + u[Comp::v()][i][j]  [k]);  // v @ jmin
    real vpf = 0.5 * (u[Comp::v()][i][j+1][k-1] + u[Comp::v()][i][j+1][k]);  // v @ jmax

//    wm = lim.limit(-vmf, u[Comp::w()][i][j+1][k],u[Comp::w()][i][j][k],u[Comp::w()][i][j-1][k]); 
//    wp = lim.limit(+vpf, u[Comp::w()][i][j-1][k],u[Comp::w()][i][j][k],u[Comp::w()][i][j+1][k]);
    real wm = 0.5 * (u[Comp::w()][i][j][k]+u[Comp::w()][i][j-1][k]);
    real wp = 0.5 * (u[Comp::w()][i][j][k]+u[Comp::w()][i][j+1][k]);

    (*conv)[Comp::w()][i][j]  [k] += a_s * vmf*wm - a_n * vpf*wp;
    }

  /* God help me, pleae */
  for_mijk(Comp::u(),i,j,k) {
    real rho = fluid()->rho(Comp::u(),i,j,k);
    (*conv)[Comp::u()][i][j][k] *= rho;
  }
  for_mijk(Comp::v(),i,j,k) {
    real rho = fluid()->rho(Comp::v(),i,j,k);
    (*conv)[Comp::v()][i][j][k] *= rho;
  }
  for_mijk(Comp::w(),i,j,k) {
    real rho = fluid()->rho(Comp::w(),i,j,k);
    (*conv)[Comp::w()][i][j][k] *= rho;
  }

  if(time->current_step() % 4000 == 0)
    boil::plot->plot(*conv, "conv", time->current_step());
}
