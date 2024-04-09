#include "momentum.h"
#include "../../../Plot/plot.h"

#define NEW_CONSERVATIVE_FORM true

/******************************************************************************/
void Momentum::convection(Vector * conv) {
/*------------+
|       3  n  |
|  f += - H   |
|       2     |
+------------*/

  u.exchange();

  for_m(m)
    for(int i=0; i<u.ni(m); i++)
      for(int j=0; j<u.nj(m); j++)
        for(int k=0; k<u.nk(m); k++) {
          (*conv)[m][i][j][k] = 0.0;
            fnew [m][i][j][k] = 0.0;
        }

  for_m(m) {
    assert( fnew.si(m)==si(m) );  assert( fnew.ei(m)==ei(m) );
    assert( fnew.sj(m)==sj(m) );  assert( fnew.ej(m)==ej(m) );
    assert( fnew.sk(m)==sk(m) );  assert( fnew.ek(m)==ek(m) );

    assert( (*conv).si(m)==si(m) );  assert( (*conv).ei(m)==ei(m) );
    assert( (*conv).sj(m)==sj(m) );  assert( (*conv).ej(m)==ej(m) );
    assert( (*conv).sk(m)==sk(m) );  assert( (*conv).ek(m)==ek(m) );
  }

  bool imin_kn, imax_kn, jmin_kn, jmax_kn, kmin_kn, kmax_kn;

  /////////
  //     //
  //  u  //
  //     //
  /////////
  Comp m = Comp::u();

  const bool imin_cont = u.bc(m).type( Dir::imin(), BndType::periodic() ) || 
                         dom->coord(m) != 0; 
  const bool imax_cont = u.bc(m).type( Dir::imax(), BndType::periodic() ) || 
                         dom->coord(m) != dom->dim(m)-1;

  imin_kn = u.bc(m).type_here( Dir::imin(), BndType::inlet() )  ||
            u.bc(m).type_here( Dir::imin(), BndType::insert() ) ||
            u.bc(m).type_here( Dir::imin(), BndType::wall() ) ;
  imax_kn = u.bc(m).type_here( Dir::imax(), BndType::inlet() )  ||
            u.bc(m).type_here( Dir::imax(), BndType::insert() ) ||
            u.bc(m).type_here( Dir::imax(), BndType::wall() ) ;
  jmin_kn = u.bc(m).type_here( Dir::jmin(), BndType::inlet() )  ||
            u.bc(m).type_here( Dir::jmin(), BndType::insert() ) ||
            u.bc(m).type_here( Dir::jmin(), BndType::wall() ) ;
  jmax_kn = u.bc(m).type_here( Dir::jmax(), BndType::inlet() )  ||
            u.bc(m).type_here( Dir::jmax(), BndType::insert() ) ||
            u.bc(m).type_here( Dir::jmax(), BndType::wall() ) ;
  kmin_kn = u.bc(m).type_here( Dir::kmin(), BndType::inlet() )  ||
            u.bc(m).type_here( Dir::kmin(), BndType::insert() ) ||
            u.bc(m).type_here( Dir::kmin(), BndType::wall() ) ;
  kmax_kn = u.bc(m).type_here( Dir::kmax(), BndType::inlet() )  ||
            u.bc(m).type_here( Dir::kmax(), BndType::insert() ) ||
            u.bc(m).type_here( Dir::kmax(), BndType::wall() ) ;

  /////////////
  //// u u ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::u(); /* transport by u */

    real a_w = dSx(m,i,j,k); 
    real a_e = dSx(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_w *= dom->ibody().fSw(m,i,j,k);
      a_e *= dom->ibody().fSe(m,i,j,k);
    } 

    real umf = 0.5 * (u[d][i-1][j][k] + u[d][i][j][k]);      // u @ imin
    real upf = 0.5 * (u[d][i+1][j][k] + u[d][i][j][k]);      // u @ imax
    
    real um = lim.limit(-umf, u[m][i+1][j][k], u[m][i][j][k], u[m][i-1][j][k]);
    real up = lim.limit(+upf, u[m][i-1][j][k], u[m][i][j][k], u[m][i+1][j][k]);

    if(i==si(m) && imin_kn) um = u[m][i-1][j][k];
    if(i==ei(m) && imax_kn) up = u[m][i+1][j][k];

    if( imin_cont && i==si(m) ) um = 0.0;
    if( imax_cont && i==ei(m) ) up = 0.0;

    (*conv)[m][i]  [j][k] += (a_w * umf*um - a_e * upf*up);
    (*conv)[m][i+1][j][k] +=  a_e * upf*up;
    (*conv)[m][i-1][j][k] -=  a_w * umf*um;
    } 

  //// exchange buffers ////
  for_mjk(m,j,k) 
   {fnew[m][ei(m)-1][j][k] = (*conv)[m][ei(m)][j][k];
    fnew[m][si(m)+1][j][k] = (*conv)[m][si(m)][j][k];}
  fnew.exchange(m,0); // in "i" direction
  for_mjk(m,j,k) 
   {(*conv)[m][ei(m)][j][k] += fnew[m][ei(m)+1][j][k];
    (*conv)[m][si(m)][j][k] += fnew[m][si(m)-1][j][k];}

  /////////////
  //// v u ////
  /////////////
  for_mijk(m,i,j,k) 
    { /* v */
    const Comp d = Comp::v(); /* transport by v */

    real a_s = dSy(m,i,j,k); 
    real a_n = dSy(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_s *= dom->ibody().fSs(m,i,j,k);
      a_n *= dom->ibody().fSn(m,i,j,k);
    } 

    real vmf = 0.5 * (u[d][i-1][j]  [k] + u[d][i][j]  [k]);  // v @ jmin
    real vpf = 0.5 * (u[d][i-1][j+1][k] + u[d][i][j+1][k]);  // v @ jmax

    real um = lim.limit(-vmf, u[m][i][j+1][k],u[m][i][j][k],u[m][i][j-1][k]);
    real up = lim.limit(+vpf, u[m][i][j-1][k],u[m][i][j][k],u[m][i][j+1][k]);
    
    if(j==sj(m) && jmin_kn) um = u[m][i][j-1][k];
    if(j==ej(m) && jmax_kn) up = u[m][i][j+1][k];

    (*conv)[m][i][j]  [k] += (a_s * vmf*um - a_n * vpf*up);
    (*conv)[m][i][j+1][k] +=  a_n * vpf*up;
    (*conv)[m][i][j-1][k] -=  a_s * vmf*um;
    } 

  //// exchange ////
  for_mik(m,i,k) 
   {fnew[m][i][ej(m)][k] = (*conv)[m][i][ej(m)+1][k];
    fnew[m][i][sj(m)][k] = (*conv)[m][i][sj(m)-1][k];}
  fnew.exchange(m,1); // in "j" direction
  for_mik(m,i,k) 
   {(*conv)[m][i][ej(m)][k] += fnew[m][i][ej(m)+1][k];
    (*conv)[m][i][sj(m)][k] += fnew[m][i][sj(m)-1][k];}

  /////////////
  //// w u ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::w(); /* transport by w */

    real a_b = dSz(m,i,j,k); 
    real a_t = dSz(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_b *= dom->ibody().fSb(m,i,j,k);
      a_t *= dom->ibody().fSt(m,i,j,k);
    } 

    real wmf = 0.5 * (u[d][i-1][j][k]   + u[d][i][j][k]);    // w @ kmin
    real wpf = 0.5 * (u[d][i-1][j][k+1] + u[d][i][j][k+1]);  // w @ kmax

    real um = lim.limit(-wmf, u[m][i][j][k+1],u[m][i][j][k],u[m][i][j][k-1]);
    real up = lim.limit(+wpf, u[m][i][j][k-1],u[m][i][j][k],u[m][i][j][k+1]);

    if(k==sk(m) && kmin_kn) um = u[m][i][j][k-1];
    if(k==ek(m) && kmax_kn) up = u[m][i][j][k+1];

    (*conv)[m][i][j][k]   += (a_b * wmf*um - a_t * wpf*up);
    (*conv)[m][i][j][k+1] +=  a_t * wpf*up;
    (*conv)[m][i][j][k-1] -=  a_b * wmf*um;
    }
  
  //// exchange ////
  for_mij(m,i,j) 
   {fnew[m][i][j][ek(m)] = (*conv)[m][i][j][ek(m)+1];
    fnew[m][i][j][sk(m)] = (*conv)[m][i][j][sk(m)-1];}
  fnew.exchange(m,2); // in "k" direction
  for_mij(m,i,j) 
   {(*conv)[m][i][j][ek(m)] += fnew[m][i][j][ek(m)+1];
    (*conv)[m][i][j][sk(m)] += fnew[m][i][j][sk(m)-1];}

  /////////
  //     //
  //  v  //
  //     //
  /////////
  m = Comp::v();

  const bool jmin_cont = u.bc(m).type( Dir::jmin(), BndType::periodic() ) || 
                         dom->coord(m) != 0; 
  const bool jmax_cont = u.bc(m).type( Dir::jmax(), BndType::periodic() ) || 
                         dom->coord(m) != dom->dim(m)-1;

  /////////////
  //// v v ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::v(); /* transport by v */

    real a_s = dSy(m,i,j,k); 
    real a_n = dSy(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_s *= dom->ibody().fSs(m,i,j,k);
      a_n *= dom->ibody().fSn(m,i,j,k);
    } 

    real vmf = 0.5 * (u[d][i][j-1][k] + u[d][i][j][k]);
    real vpf = 0.5 * (u[d][i][j+1][k] + u[d][i][j][k]);
    
    real vm = lim.limit(-vmf, u[m][i][j+1][k],u[m][i][j][k],u[m][i][j-1][k]);
    real vp = lim.limit(+vpf, u[m][i][j-1][k],u[m][i][j][k],u[m][i][j+1][k]);

    if(j==sj(m) && jmin_kn) vm = u[m][i][j-1][k];
    if(j==ej(m) && jmax_kn) vp = u[m][i][j+1][k];

    if( jmin_cont && j==sj(m) ) vm = 0.0;
    if( jmax_cont && j==ej(m) ) vp = 0.0;

    (*conv)[m][i][j]  [k] += (a_s * vmf*vm - a_n * vpf*vp);
    (*conv)[m][i][j+1][k] +=  a_n * vpf*vp;
    (*conv)[m][i][j-1][k] -=  a_s * vmf*vm;
    }

  //// exchange buffers ////
  for_mik(m,i,k) 
   {fnew[m][i][ej(m)-1][k] = (*conv)[m][i][ej(m)][k];
    fnew[m][i][sj(m)+1][k] = (*conv)[m][i][sj(m)][k];}
  fnew.exchange(m,1); // in "j" direction
  for_mik(m,i,k) 
   {(*conv)[m][i][ej(m)][k] += fnew[m][i][ej(m)+1][k];
    (*conv)[m][i][sj(m)][k] += fnew[m][i][sj(m)-1][k];}

  /////////////
  //// u v ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::u(); /* transport by u */

    real a_w = dSx(m,i,j,k); 
    real a_e = dSx(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_w *= dom->ibody().fSw(m,i,j,k);
      a_e *= dom->ibody().fSe(m,i,j,k);
    } 

    real umf = 0.5 * (u[d][i]  [j][k] + u[d][i]  [j-1][k]);  // u @ imin
    real upf = 0.5 * (u[d][i+1][j][k] + u[d][i+1][j-1][k]);  // u @ imax

    real vm = lim.limit(-umf, u[m][i+1][j][k],u[m][i][j][k],u[m][i-1][j][k]);
    real vp = lim.limit(+upf, u[m][i-1][j][k],u[m][i][j][k],u[m][i+1][j][k]);

    if(i==si(m) && imin_kn) vm = u[m][i-1][j][k];
    if(i==ei(m) && imax_kn) vp = u[m][i+1][j][k];

    (*conv)[m][i]  [j][k] += (a_w * umf*vm - a_e * upf*vp);
    (*conv)[m][i+1][j][k] +=  a_e * upf*vp;
    (*conv)[m][i-1][j][k] -=  a_w * umf*vm;
    }

  //// exchange ////
  for_mjk(m,j,k) 
   {fnew[m][ei(m)][j][k] = (*conv)[m][ei(m)+1][j][k];
    fnew[m][si(m)][j][k] = (*conv)[m][si(m)-1][j][k];}
  fnew.exchange(m,0); // in "i" direction
  for_mjk(m,j,k) 
   {(*conv)[m][ei(m)][j][k] += fnew[m][ei(m)+1][j][k];
    (*conv)[m][si(m)][j][k] += fnew[m][si(m)-1][j][k];}

  /////////////
  //// w v ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::w(); /* transport by w */

    real a_b = dSz(m,i,j,k); 
    real a_t = dSz(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_b *= dom->ibody().fSb(m,i,j,k);
      a_t *= dom->ibody().fSt(m,i,j,k);
    } 

    real wmf = 0.5 * (u[d][i][j][k]   + u[d][i][j-1][k]);    // w @ kmin
    real wpf = 0.5 * (u[d][i][j][k+1] + u[d][i][j-1][k+1]);  // w @ kmax

    real vm = lim.limit(-wmf, u[m][i][j][k+1],u[m][i][j][k],u[m][i][j][k-1]);
    real vp = lim.limit(+wpf, u[m][i][j][k-1],u[m][i][j][k],u[m][i][j][k+1]);

    if(k==sk(m) && kmin_kn) vm = u[m][i][j][k-1];
    if(k==ek(m) && kmax_kn) vp = u[m][i][j][k+1];

    (*conv)[m][i][j][k]   += (a_b * wmf*vm - a_t * wpf*vp);
    (*conv)[m][i][j][k+1] +=  a_t * wpf*vp;
    (*conv)[m][i][j][k-1] -=  a_b * wmf*vm;
    }

  //// exchange ////
  for_mij(m,i,j) 
   {fnew[m][i][j][ek(m)] = (*conv)[m][i][j][ek(m)+1];
    fnew[m][i][j][sk(m)] = (*conv)[m][i][j][sk(m)-1];}
  fnew.exchange(m,2); // in "k" direction
  for_mij(m,i,j) 
   {(*conv)[m][i][j][ek(m)] += fnew[m][i][j][ek(m)+1];
    (*conv)[m][i][j][sk(m)] += fnew[m][i][j][sk(m)-1];}
  
  /////////
  //     //
  //  w  //
  //     //
  /////////
  m = Comp::w();

  const bool kmin_cont = u.bc(m).type( Dir::kmin(), BndType::periodic() ) || 
                         dom->coord(m) != 0; 
  const bool kmax_cont = u.bc(m).type( Dir::kmax(), BndType::periodic() ) || 
                         dom->coord(m) != dom->dim(m)-1;

  /////////////
  //// w w ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::w(); /* transport by w */

    real a_b = dSz(m,i,j,k); 
    real a_t = dSz(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_b *= dom->ibody().fSb(m,i,j,k);
      a_t *= dom->ibody().fSt(m,i,j,k);
    } 

    real wmf = 0.5 * (u[d][i][j][k-1] + u[d][i][j][k]);      // w @ kmin
    real wpf = 0.5 * (u[d][i][j][k+1] + u[d][i][j][k]);      // w @ kmax
    
    real wm = lim.limit(-wmf, u[m][i][j][k+1],u[m][i][j][k],u[m][i][j][k-1]);
    real wp = lim.limit(+wpf, u[m][i][j][k-1],u[m][i][j][k],u[m][i][j][k+1]);

    if(k==sk(m) && kmin_kn) wm = u[m][i][j][k-1];
    if(k==ek(m) && kmax_kn) wp = u[m][i][j][k+1];

    if( kmin_cont && k==sk(m) ) wm = 0.0;
    if( kmax_cont && k==ek(m) ) wp = 0.0;

    (*conv)[m][i][j][k]   += (a_b * wmf*wm - a_t * wpf*wp);
    (*conv)[m][i][j][k+1] +=  a_t * wpf*wp;
    (*conv)[m][i][j][k-1] -=  a_b * wmf*wm;
    }

  //// exchange buffers ////
  for_mij(m,i,j) 
   {fnew[m][i][j][ek(m)-1] = (*conv)[m][i][j][ek(m)];
    fnew[m][i][j][sk(m)+1] = (*conv)[m][i][j][sk(m)];}
  fnew.exchange(m,2); // in "k" direction
  for_mij(m,i,j) 
   {(*conv)[m][i][j][ek(m)] += fnew[m][i][j][ek(m)+1];
    (*conv)[m][i][j][sk(m)] += fnew[m][i][j][sk(m)-1];}

  /////////////
  //// u w ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::u(); /* transport by u */

    real a_w = dSx(m,i,j,k); 
    real a_e = dSx(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_w *= dom->ibody().fSw(m,i,j,k);
      a_e *= dom->ibody().fSe(m,i,j,k);
    } 

    real umf = 0.5 * (u[d][i]  [j][k-1] + u[d][i]  [j][k]);  // u @ imin
    real upf = 0.5 * (u[d][i+1][j][k-1] + u[d][i+1][j][k]);  // u @ imax

    real wm = lim.limit(-umf, u[m][i+1][j][k],u[m][i][j][k],u[m][i-1][j][k]);
    real wp = lim.limit(+upf, u[m][i-1][j][k],u[m][i][j][k],u[m][i+1][j][k]);

    if(i==si(m) && imin_kn) wm = u[m][i-1][j][k];
    if(i==ei(m) && imax_kn) wp = u[m][i+1][j][k];

    (*conv)[m][i]  [j][k] += (a_w * umf*wm - a_e * upf*wp);
    (*conv)[m][i+1][j][k] +=  a_e * upf*wp;
    (*conv)[m][i-1][j][k] -=  a_w * umf*wm;
    }

  //// exchange ////
  for_mjk(m,j,k) 
   {fnew[m][ei(m)][j][k] = (*conv)[m][ei(m)+1][j][k];
    fnew[m][si(m)][j][k] = (*conv)[m][si(m)-1][j][k];}
  fnew.exchange(m,0); // in "i" direction
  for_mjk(m,j,k) 
   {(*conv)[m][ei(m)][j][k] += fnew[m][ei(m)+1][j][k];
    (*conv)[m][si(m)][j][k] += fnew[m][si(m)-1][j][k];}

  /////////////
  //// v w ////
  /////////////
  for_mijk(m,i,j,k) 
    { 
    const Comp d = Comp::v(); /* transport by v */

    real a_s = dSy(m,i,j,k); 
    real a_n = dSy(m,i,j,k); 
    if( dom->ibody().cut(m,i,j,k) ) {
      a_s *= dom->ibody().fSs(m,i,j,k);
      a_n *= dom->ibody().fSn(m,i,j,k);
    } 

    real vmf = 0.5 * (u[d][i][j]  [k-1] + u[d][i][j]  [k]);  // v @ jmin
    real vpf = 0.5 * (u[d][i][j+1][k-1] + u[d][i][j+1][k]);  // v @ jmax

    real wm = lim.limit(-vmf, u[m][i][j+1][k],u[m][i][j][k],u[m][i][j-1][k]);
    real wp = lim.limit(+vpf, u[m][i][j-1][k],u[m][i][j][k],u[m][i][j+1][k]);

    if(j==sj(m) && jmin_kn) wm = u[m][i][j-1][k];
    if(j==ej(m) && jmax_kn) wp = u[m][i][j+1][k];

    (*conv)[m][i][j]  [k] += (a_s * vmf*wm - a_n * vpf*wp);
    (*conv)[m][i][j+1][k] +=  a_n * vpf*wp;
    (*conv)[m][i][j-1][k] -=  a_s * vmf*wm;
    }

  //// exchange ////
  for_mik(m,i,k) 
   {fnew[m][i][ej(m)][k] = (*conv)[m][i][ej(m)+1][k];
    fnew[m][i][sj(m)][k] = (*conv)[m][i][sj(m)-1][k];}
  fnew.exchange(m,1); // in "j" direction
  for_mik(m,i,k) 
   {(*conv)[m][i][ej(m)][k] += fnew[m][i][ej(m)+1][k];
    (*conv)[m][i][sj(m)][k] += fnew[m][i][sj(m)-1][k];}

  /* yohei's correction in new form */
  /* beware: it doesn't take care of immersed boundaries */
  #if NEW_CONSERVATIVE_FORM
  /* u.div(uvw) */
  m=Comp::u();
  for_mijk(m,i,j,k) {
    real divuc= - u.domain()->dSx(i  ,j,k)*u[Comp::u()][i]  [j]  [k]
                + u.domain()->dSx(i  ,j,k)*u[Comp::u()][i+1][j]  [k]
                - u.domain()->dSy(i  ,j,k)*u[Comp::v()][i]  [j]  [k]
                + u.domain()->dSy(i  ,j,k)*u[Comp::v()][i]  [j+1][k]
                - u.domain()->dSz(i  ,j,k)*u[Comp::w()][i]  [j]  [k]
                + u.domain()->dSz(i  ,j,k)*u[Comp::w()][i]  [j]  [k+1];
    real divum= - u.domain()->dSx(i-1,j,k)*u[Comp::u()][i-1][j]  [k]
                + u.domain()->dSx(i-1,j,k)*u[Comp::u()][i]  [j]  [k]
                - u.domain()->dSy(i-1,j,k)*u[Comp::v()][i-1][j]  [k]
                + u.domain()->dSy(i-1,j,k)*u[Comp::v()][i-1][j+1][k]
                - u.domain()->dSz(i-1,j,k)*u[Comp::w()][i-1][j]  [k]
                + u.domain()->dSz(i-1,j,k)*u[Comp::w()][i-1][j]  [k+1];
    (*conv)[m][i][j][k] += u[m][i][j][k] * 0.5 * (divuc+divum);
  }
  /* v.div(uvw) */
  m=Comp::v();
  for_mijk(m,i,j,k) {
    real divuc= - u.domain()->dSx(i,j  ,k)*u[Comp::u()][i]  [j]  [k]
                + u.domain()->dSx(i,j  ,k)*u[Comp::u()][i+1][j]  [k]
                - u.domain()->dSy(i,j  ,k)*u[Comp::v()][i]  [j]  [k]
                + u.domain()->dSy(i,j  ,k)*u[Comp::v()][i]  [j+1][k]
                - u.domain()->dSz(i,j  ,k)*u[Comp::w()][i]  [j]  [k]
                + u.domain()->dSz(i,j  ,k)*u[Comp::w()][i]  [j]  [k+1];
    real divum= - u.domain()->dSx(i,j-1,k)*u[Comp::u()][i]  [j-1][k]
                + u.domain()->dSx(i,j-1,k)*u[Comp::u()][i+1][j-1][k]
                - u.domain()->dSy(i,j-1,k)*u[Comp::v()][i]  [j-1][k]
                + u.domain()->dSy(i,j-1,k)*u[Comp::v()][i]  [j  ][k]
                - u.domain()->dSz(i,j-1,k)*u[Comp::w()][i]  [j-1][k]
                + u.domain()->dSz(i,j-1,k)*u[Comp::w()][i]  [j-1][k+1];
    (*conv)[m][i][j][k] += u[m][i][j][k] * 0.5 * (divuc+divum);
  }
  /* w.div(uvw) */
  m=Comp::w();
  for_mijk(m,i,j,k) {
    real divuc= - u.domain()->dSx(i,j,k  )*u[Comp::u()][i]  [j]  [k]
                + u.domain()->dSx(i,j,k  )*u[Comp::u()][i+1][j]  [k]
                - u.domain()->dSy(i,j,k  )*u[Comp::v()][i]  [j]  [k]
                + u.domain()->dSy(i,j,k  )*u[Comp::v()][i]  [j+1][k]
                - u.domain()->dSz(i,j,k  )*u[Comp::w()][i]  [j]  [k]
                + u.domain()->dSz(i,j,k  )*u[Comp::w()][i]  [j]  [k+1];
    real divum= - u.domain()->dSx(i,j,k-1)*u[Comp::u()][i]  [j]  [k-1]
                + u.domain()->dSx(i,j,k-1)*u[Comp::u()][i+1][j]  [k-1]
                - u.domain()->dSy(i,j,k-1)*u[Comp::v()][i]  [j]  [k-1]
                + u.domain()->dSy(i,j,k-1)*u[Comp::v()][i]  [j+1][k-1]
                - u.domain()->dSz(i,j,k-1)*u[Comp::w()][i]  [j]  [k-1]
                + u.domain()->dSz(i,j,k-1)*u[Comp::w()][i]  [j]  [k];
    (*conv)[m][i][j][k] += u[m][i][j][k] * 0.5 * (divuc+divum);
  }
  #endif

  /* God help me, pleae */
  for_m(m)
    for_mijk(m,i,j,k) (*conv)[m][i][j][k] *= fluid()->rho(m,i,j,k);
  
  //if(time->current_step() % 1000 == 0)
  //  boil::plot->plot(*conv, "conv", time->current_step());
}
