#include "momentum.h"

/******************************************************************************/
void Momentum::create_system(const Scalar * mu_eddy) {
/*-----------------------+ 
|  create system matrix  | -> should BndType::outlet() be checked at all?
+-----------------------*/

  real mu, rho;
  real mu_w, mu_e, mu_s, mu_n, mu_b, mu_t;
  real a_x, a_y, a_z;

  /*------------------------------+
  |  coefficient due to innertia  |
  +------------------------------*/
  Comp m = Comp::u();
  for_mijk(m,i,j,k) {
    /* linear */
    rho = 0.5 * (fluid()->rho(i,j,k) + fluid()->rho(i-1,j,k));
    A[~m]->c[i][j][k] = dV(m,i,j,k) * rho * time->dti();
  }
  m = Comp::v();
  for_mijk(m,i,j,k) {
    /* linear */
    rho = 0.5 * (fluid()->rho(i,j,k) + fluid()->rho(i,j-1,k));
    A[~m]->c[i][j][k] = dV(m,i,j,k) * rho * time->dti();
  }
  m = Comp::w();
  for_mijk(m,i,j,k) {
    /* linear */
    rho = 0.5 * (fluid()->rho(i,j,k) + fluid()->rho(i,j,k-1));
    A[~m]->c[i][j][k] = dV(m,i,j,k) * rho * time->dti();
  }
 
  /* correct for volumes in immersed boundary */
  for_m(m) {
    if(dom->ibody().nccells(m) > 0) {
      for(int cc=0; cc<dom->ibody().nccells(m); cc++)
        if( dom->ibody().cut(m,cc) ) {
          int i,j,k;
          dom->ibody().ijk(m,cc,&i,&j,&k); 
          A[~m]->c[i][j][k] *= dom->ibody().fV(m,cc); 
        }
    }
  }

  /*--------------------------------+
  |  coefficients due to viscosity  |
  +--------------------------------*/

  /* get time stepping coefficient */
  real tsc = diff_ts.N();
  assert( tsc > 0.0 );

  /*-------------------------------------------------------------+
  |  Following is crude and dangerous code for immersed boudary. |
  |  mu inside solid is directly used.                           |
  |  If color function is extrapolated from fliud to solid,      |
  |  then the code returns correct matrix.                       |
  +-------------------------------------------------------------*/

  /////////
  //     //
  //  u  //
  //     //
  /////////
  m = Comp::u();
  for_mijk(m,i,j,k) {
    /*----------------------------------------+ 
    |  coefficients in i direction (w and e)  |
    +----------------------------------------*/
    a_x  = dSx(m,i,j,k); // dyc(m,j) * dzc(m,k);
    if( !mu_eddy )
     {mu_w = fluid()->mu(i-1,j,k);
      mu_e = fluid()->mu(i  ,j,k);}
    else
     {mu_w = fluid()->mu(i-1,j,k) + (*mu_eddy)[i-1][j][k];
      mu_e = fluid()->mu(i  ,j,k) + (*mu_eddy)[i]  [j][k];}
    A[~m]->w[i][j][k] = tsc * mu_w * a_x / dxw(m,i);
    A[~m]->e[i][j][k] = tsc * mu_e * a_x / dxe(m,i);
    
    /*----------------------------------------+ 
    |  coefficients in j direction (s and n)  |
    +----------------------------------------*/
    a_y  = dSy(m,i,j,k); //dxc(m,i) * dzc(m,k);
    if( !mu_eddy )
     {mu_s = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j-1,k) 
                   + fluid()->mu(i  ,j-1,k));
      mu_n = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j+1,k) 
                   + fluid()->mu(i  ,j+1,k));}
    else
     {mu_s = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j-1,k) + (*mu_eddy)[i-1][j-1][k]
                   + fluid()->mu(i  ,j-1,k) + (*mu_eddy)[i]  [j-1][k]);
      mu_n = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j+1,k) + (*mu_eddy)[i-1][j+1][k]
                   + fluid()->mu(i  ,j+1,k) + (*mu_eddy)[i]  [j+1][k]);}
    A[~m]->s[i][j][k] = tsc * mu_s * a_y / dys(m,j);
    A[~m]->n[i][j][k] = tsc * mu_n * a_y / dyn(m,j);
    
    /*----------------------------------------+ 
    |  coefficients in k direction (b and t)  |
    +----------------------------------------*/
    a_z  = dSz(m,i,j,k); //dyc(m,j) * dxc(m,i);
    if( !mu_eddy )
     {mu_b = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j,k-1) 
                   + fluid()->mu(i  ,j,k-1));
      mu_t = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j,k+1) 
                   + fluid()->mu(i  ,j,k+1));}
    else
     {mu_b = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j,k-1) + (*mu_eddy)[i-1][j][k-1]
                   + fluid()->mu(i  ,j,k-1) + (*mu_eddy)[i]  [j][k-1]);
      mu_t = 0.25 * (mu_w + mu_e)
           + 0.25 * (fluid()->mu(i-1,j,k+1) + (*mu_eddy)[i-1][j][k+1]
                   + fluid()->mu(i  ,j,k+1) + (*mu_eddy)[i]  [j][k+1]);}
    A[~m]->b[i][j][k] = tsc * mu_b * a_z / dzb(m,k);
    A[~m]->t[i][j][k] = tsc * mu_t * a_z / dzt(m,k);
    
  }

  /////////
  //     //
  //  v  //
  //     //
  /////////
  m = Comp::v();
  for_mijk(m,i,j,k) {
    /*----------------------------------------+ 
    |  coefficients in j direction (s and n)  |
    +----------------------------------------*/
    a_y  = dSy(m,i,j,k); //dxc(m,i) * dzc(m,k);
    if( !mu_eddy )
     {mu_s = fluid()->mu(i,j-1,k);
      mu_n = fluid()->mu(i,j  ,k);}
    else
     {mu_s = fluid()->mu(i,j-1,k) + (*mu_eddy)[i][j-1][k];
      mu_n = fluid()->mu(i,j  ,k) + (*mu_eddy)[i][j]  [k];}
    A[~m]->s[i][j][k] = tsc * mu_s * a_y / dys(m,j);
    A[~m]->n[i][j][k] = tsc * mu_n * a_y / dyn(m,j);
    
    /*----------------------------------------+ 
    |  coefficients in i direction (w and e)  |
    +----------------------------------------*/
    a_x  = dSx(m,i,j,k); //dyc(m,j) * dzc(m,k);
    if( !mu_eddy )
     {mu_w = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i-1,j-1,k) 
                   + fluid()->mu(i-1,j  ,k));
      mu_e = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i+1,j-1,k) 
                   + fluid()->mu(i+1,j  ,k));}
    else
     {mu_w = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i-1,j-1,k) + (*mu_eddy)[i-1][j-1][k]
                   + fluid()->mu(i-1,j  ,k) + (*mu_eddy)[i-1][j]  [k]);
      mu_e = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i+1,j-1,k) + (*mu_eddy)[i+1][j-1][k] 
                   + fluid()->mu(i+1,j  ,k) + (*mu_eddy)[i+1][j]  [k]);}
    A[~m]->w[i][j][k] = tsc * mu_w * a_x / dxw(m,i);
    A[~m]->e[i][j][k] = tsc * mu_e * a_x / dxe(m,i);
    
    /*----------------------------------------+ 
    |  coefficients in k direction (b and t)  |
    +----------------------------------------*/
    a_z  = dSz(m,i,j,k); //dyc(m,j) * dxc(m,i);
    if( !mu_eddy )
     {mu_b = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i,j-1,k-1) 
                   + fluid()->mu(i,j  ,k-1));
      mu_t = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i,j-1,k+1) 
                   + fluid()->mu(i,j  ,k+1));}
    else
     {mu_b = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i,j-1,k-1) + (*mu_eddy)[i][j-1][k-1] 
                   + fluid()->mu(i,j  ,k-1) + (*mu_eddy)[i][j]  [k-1]);
      mu_t = 0.25 * (mu_s + mu_n)
           + 0.25 * (fluid()->mu(i,j-1,k+1) + (*mu_eddy)[i][j-1][k+1] 
                   + fluid()->mu(i,j  ,k+1) + (*mu_eddy)[i][j]  [k+1]);}
    A[~m]->b[i][j][k] = tsc * mu_b * a_z / dzb(m,k);
    A[~m]->t[i][j][k] = tsc * mu_t * a_z / dzt(m,k);
    
  }

  /////////
  //     //
  //  w  //
  //     //
  /////////
  m = Comp::w();
  for_mijk(m,i,j,k) {
    /*----------------------------------------+ 
    |  coefficients in k direction (b and t)  |
    +----------------------------------------*/
    a_z  = dSz(m,i,j,k); //dyc(m,j) * dxc(m,i);
    if( !mu_eddy )
     {mu_b = fluid()->mu(i,j,k-1);
      mu_t = fluid()->mu(i,j,k  );}
    else
     {mu_b = fluid()->mu(i,j,k-1) + (*mu_eddy)[i][j][k-1];
      mu_t = fluid()->mu(i,j,k  ) + (*mu_eddy)[i][j][k]  ;}
    A[~m]->b[i][j][k] = tsc * mu_b * a_z / dzb(m,k);
    A[~m]->t[i][j][k] = tsc * mu_t * a_z / dzt(m,k);
    
    /*----------------------------------------+ 
    |  coefficients in i direction (w and e)  |
    +----------------------------------------*/
    a_x  = dSx(m,i,j,k); //dyc(m,j) * dzc(m,k);
    if( !mu_eddy )
     {mu_w = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i-1,j,k-1) 
                   + fluid()->mu(i-1,j,k  ));
      mu_e = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i+1,j,k-1) 
                   + fluid()->mu(i+1,j,k  ));}
    else
     {mu_w = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i-1,j,k-1)+ (*mu_eddy)[i-1][j][k-1] 
                   + fluid()->mu(i-1,j,k  )+ (*mu_eddy)[i-1][j][k]  );
      mu_e = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i+1,j,k-1)+ (*mu_eddy)[i+1][j][k-1] 
                   + fluid()->mu(i+1,j,k  )+ (*mu_eddy)[i+1][j][k]  );}
    A[~m]->w[i][j][k] = tsc * mu_w * a_x / dxw(m,i);
    A[~m]->e[i][j][k] = tsc * mu_e * a_x / dxe(m,i);
    
    /*----------------------------------------+ 
    |  coefficients in j direction (s and n)  |
    +----------------------------------------*/
    a_y  = dSy(m,i,j,k); //dxc(m,i) * dzc(m,k);
    if( !mu_eddy )
     {mu_s = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i,j-1,k-1) 
                   + fluid()->mu(i,j-1,k  ));
      mu_n = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i,j+1,k-1) 
                   + fluid()->mu(i,j+1,k  ));}
    else
     {mu_s = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i,j-1,k-1)+ (*mu_eddy)[i][j-1][k-1] 
                   + fluid()->mu(i,j-1,k  )+ (*mu_eddy)[i][j-1][k]  );
      mu_n = 0.25 * (mu_b + mu_t)
           + 0.25 * (fluid()->mu(i,j+1,k-1)+ (*mu_eddy)[i][j+1][k-1] 
                   + fluid()->mu(i,j+1,k  )+ (*mu_eddy)[i][j+1][k]  );}
    A[~m]->s[i][j][k] = tsc * mu_s * a_y / dys(m,j);
    A[~m]->n[i][j][k] = tsc * mu_n * a_y / dyn(m,j);
    
  }

  /*-------------------------------+
  |  a "touch" from immersed body  |
  +-------------------------------*/
  for_m(m) 
    if(dom->ibody().nccells(m) > 0) {
    
      for(int cc=0; cc<dom->ibody().nccells(m); cc++) {
        int i,j,k;
        dom->ibody().ijk(m,cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);


        /* define some aliases */
        const real fSw = dom->ibody().fSw(m,cc);
        const real fSe = dom->ibody().fSe(m,cc);
        const real fSs = dom->ibody().fSs(m,cc);
        const real fSn = dom->ibody().fSn(m,cc);
        const real fSb = dom->ibody().fSb(m,cc);
        const real fSt = dom->ibody().fSt(m,cc);

        const real fdxw = dom->ibody().fdxw(m,cc);
        const real fdxe = dom->ibody().fdxe(m,cc);
        const real fdys = dom->ibody().fdys(m,cc);
        const real fdyn = dom->ibody().fdyn(m,cc);
        const real fdzb = dom->ibody().fdzb(m,cc);
        const real fdzt = dom->ibody().fdzt(m,cc);

        /* cell i,j,k is "cut" and "on" */
        if( dom->ibody().on(m,i,j,k) ) {

          if( dom->ibody().cut(m,i-1,j,k) ) A[~m]->w[i][j][k] *= fSw;
          if( dom->ibody().off(m,i-1,j,k) ) A[~m]->w[i][j][k] /= fdxw;

          if( dom->ibody().cut(m,i+1,j,k) ) A[~m]->e[i][j][k] *= fSe;
          if( dom->ibody().off(m,i+1,j,k) ) A[~m]->e[i][j][k] /= fdxe;

          if( dom->ibody().cut(m,i,j-1,k) ) A[~m]->s[i][j][k] *= fSs;
          if( dom->ibody().off(m,i,j-1,k) ) A[~m]->s[i][j][k] /= fdys;

          if( dom->ibody().cut(m,i,j+1,k) ) A[~m]->n[i][j][k] *= fSn;
          if( dom->ibody().off(m,i,j+1,k) ) A[~m]->n[i][j][k] /= fdyn;

          if( dom->ibody().cut(m,i,j,k-1) ) A[~m]->b[i][j][k] *= fSb;
          if( dom->ibody().off(m,i,j,k-1) ) A[~m]->b[i][j][k] /= fdzb;

          if( dom->ibody().cut(m,i,j,k+1) ) A[~m]->t[i][j][k] *= fSt;
          if( dom->ibody().off(m,i,j,k+1) ) A[~m]->t[i][j][k] /= fdzt;

        /* cell i,j,k is "cut" but "off" */
        } else {
          if(fdxe != 1.0) A[~m]->w[i+1][j][k] /= (1.0-fdxe);
          if(fdxw != 1.0) A[~m]->e[i-1][j][k] /= (1.0-fdxw);
           
          if(fdys != 1.0) A[~m]->n[i][j-1][k] /= (1.0-fdys);
          if(fdyn != 1.0) A[~m]->s[i][j+1][k] /= (1.0-fdyn);

          if(fdzb != 1.0) A[~m]->t[i][j][k-1] /= (1.0-fdzb);
          if(fdzt != 1.0) A[~m]->b[i][j][k+1] /= (1.0-fdzt);
        }
      }
    }


  for_m(m) { 

    /*----------------------+ 
    |  symmetry and outlet  |
    +----------------------*/
    for( int b=0; b<u.bc(m).count(); b++ ) {

      if( u.bc(m).type(b) == BndType::symmetry() ||
          u.bc(m).type(b) == BndType::outlet() ) {

        Dir d = u.bc(m).direction(b);

        if( d == Dir::imin() ) 
          for_vjk(u.bc(m).at(b),j,k) 
            A[~m]->w[si(m)][j][k]  = 0.0;

        if( d == Dir::imax() ) 
          for_vjk(u.bc(m).at(b),j,k) 
            A[~m]->e[ei(m)][j][k]  = 0.0;

        if( d == Dir::jmin() ) 
          for_vik(u.bc(m).at(b),i,k) 
            A[~m]->s[i][sj(m)][k]  = 0.0;

        if( d == Dir::jmax() ) 
          for_vik(u.bc(m).at(b),i,k) 
            A[~m]->n[i][ej(m)][k]  = 0.0;

        if( d == Dir::kmin() ) 
          for_vij(u.bc(m).at(b),i,j) 
            A[~m]->b[i][j][sk(m)]  = 0.0;

        if( d == Dir::kmax() ) 
          for_vij(u.bc(m).at(b),i,j) 
            A[~m]->t[i][j][ek(m)]  = 0.0;

      } // symmetry or outlet
    } // b

    for_mijk(m,i,j,k) 
      A[~m]->c[i][j][k] += A[~m]->w[i][j][k] + A[~m]->e[i][j][k]
                        +  A[~m]->s[i][j][k] + A[~m]->n[i][j][k]
                        +  A[~m]->b[i][j][k] + A[~m]->t[i][j][k];

    /*-------------------------------+
    |  a "touch" from immersed body  |
    +-------------------------------*/
    if(dom->ibody().nccells() > 0) 
      for_mijk(m,i,j,k) 
        if( dom->ibody().off(m,i,j,k) ) {
          A[~m]->w[i][j][k] = 0.0;
          A[~m]->e[i][j][k] = 0.0;
          A[~m]->s[i][j][k] = 0.0;
          A[~m]->n[i][j][k] = 0.0;
          A[~m]->b[i][j][k] = 0.0;
          A[~m]->t[i][j][k] = 0.0;
          A[~m]->c[i][j][k] = 1.0;
        }

    /*------------------------------+ 
    |  compute central coefficient  |
    +------------------------------*/
    for_mijk(m,i,j,k) 
      A[~m]->ci[i][j][k] = 1.0 / A[~m]->c[i][j][k];

  } // through m, vector components

  for_m(m){
    A[~m]->c.exchange();
    A[~m]->w.exchange();
    A[~m]->e.exchange();
    A[~m]->s.exchange();
    A[~m]->n.exchange();
    A[~m]->b.exchange();
    A[~m]->t.exchange();
    A[~m]->ci.exchange();
  }

}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_create_system.cpp,v 1.34 2014/08/06 08:42:43 sato Exp $'/
+-----------------------------------------------------------------------------*/
