#include "momentum.h"

/******************************************************************************/
real Momentum::bulk(const Comp & m, const real & coord) const {

  assert( m==Comp::i() || m==Comp::j() || m==Comp::k() );

  if(m==Comp::i())      return bulk_i(coord);
  else if(m==Comp::j()) return bulk_j(coord);
  else                  return bulk_k(coord);
}


/******************************************************************************/
real Momentum::bulk_i(const real & xp) const {
/*---------------------------------------+ 
|  compute bulk velocity at specified x  |
+---------------------------------------*/

  int  ip = -1;
  real A  = 0.0;
  real b  = 0.0;

  /*---------+
  |  find i  |
  +---------*/
  for(int i=1; i<dom->ni()-2; i++) 
    if(xp > dom->xc(i)  &&  xp <= dom->xc(i+1))
     {ip = i+1;
      break;}

  /*-------------------+
  |  if plane is here  |
  +-------------------*/
  Comp m = Comp::u();

  if(ip != -1) {
    /* regular cells */
    for_mjk(m,j,k) 
      if( dom->ibody().on(ip,j,k) && !dom->ibody().cut(ip,j,k) ) {
        real a_x  = dSx(m,ip,j,k); 
        b += a_x * u[m][ip][j][k];
        A += a_x;
      }
    /* cut cells */
    for(int cc=0; cc<dom->ibody().nccells(); cc++) { /* scalar dim. ok */
      int i,j,k;
      dom->ibody().ijk(cc,&i,&j,&k);
      if(i==ip) {
        real a_x  = dSx(m,ip,j,k) * dom->ibody().fSw(cc); /* scalar dim. ok */ 
        b += a_x * u[m][ip][j][k];
        A += a_x;
      }
    }
  }

  boil::cart.sum_real(&A);
  boil::cart.sum_real(&b);

  assert(A > 0.0);

  return b/A;
}

/******************************************************************************/
real Momentum::bulk_j(const real & yp) const {
/*---------------------------------------+ 
|  compute bulk velocity at specified y  |
+---------------------------------------*/

  int  jp = -1;
  real A  = 0.0;
  real b  = 0.0;

  /*---------+
  |  find j  |
  +---------*/
  for(int j=1; j<dom->nj()-2; j++) 
    if(yp > dom->yc(j)  &&  yp <= dom->yc(j+1))
     {jp = j+1;
      break;}

  /*-------------------+
  |  if plane is here  |
  +-------------------*/
  Comp m = Comp::v();

  if(jp != -1) {
    /* regular cells */
    for_mik(m,i,k) 
      if( dom->ibody().on(i,jp,k) && !dom->ibody().cut(i,jp,k) ) {
        real a_y  = dSy(m,i,jp,k); //dxc(m,i) * dzc(m,k);
        b += a_y * u[m][i][jp][k];
        A += a_y;
      }
    /* cut cells */
    for(int cc=0; cc<dom->ibody().nccells(); cc++) { /* scalar dim. ok */
      int i,j,k;
      dom->ibody().ijk(cc,&i,&j,&k);
      if(j==jp) {
        real a_y  = dSy(m,i,jp,k) * dom->ibody().fSs(cc); /* scalar dim. ok */ 
        b += a_y * u[m][i][jp][k];
        A += a_y;
      }
    }
  }

  boil::cart.sum_real(&A);
  boil::cart.sum_real(&b);

  assert(A > 0.0);

  return b/A;
}

/******************************************************************************/
real Momentum::bulk_k(const real & zp) const {
/*---------------------------------------+ 
|  compute bulk velocity at specified z  |
+---------------------------------------*/

  int  kp = -1;
  real A  = 0.0;
  real b  = 0.0;

  /*---------+
  |  find k  |
  +---------*/
  for(int k=1; k<dom->nk()-2; k++) 
    if(zp > dom->zc(k)  &&  zp <= dom->zc(k+1))
     {kp = k+1;
      break;}

  /*-------------------+
  |  if plane is here  |
  +-------------------*/
  Comp m = Comp::w();

  if(kp != -1) {
    /* regular cells */
    for_mij(m,i,j) 
      if( dom->ibody().on(i,j,kp) && !dom->ibody().cut(i,j,kp) ) {
        real a_z  = dSz(m,i,j,kp); //dxc(m,i) * dzc(m,k);
        b += a_z * u[m][i][j][kp];
        A += a_z;
      }
    /* cut cells */
    for(int cc=0; cc<dom->ibody().nccells(); cc++) { /* scalar dim. ok */
      int i,j,k;
      dom->ibody().ijk(cc,&i,&j,&k);
      if(k==kp) {
        real a_z  = dSz(m,i,j,kp) * dom->ibody().fSb(cc); /* scalar dim. ok */ 
        b += a_z * u[m][i][j][kp];
        A += a_z;
      }
    }
  }

  boil::cart.sum_real(&A);
  boil::cart.sum_real(&b);

  assert(A > 0.0);

  return b/A;
}
