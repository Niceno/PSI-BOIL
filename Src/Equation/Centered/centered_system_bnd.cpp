#include "centered.h"

/***************************************************************************//**
*  \brief Corrects system matrix \f$ [A] \f$ at boundaries.
*******************************************************************************/
void Centered::create_system_bnd(const Property * f_prop) {

  /*--------------------------------------+
  |  initialize, that is quite important  |
  +--------------------------------------*/
  fbnd = 0.0;
  phi.bnd_update();

  /*----------------------+ 
  |  symmetry and outlet  |
  +----------------------*/
  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type(b) == BndType::neumann() ||
        phi.bc().type(b) == BndType::symmetry() ||
        phi.bc().type(b) == BndType::pseudo() ||
        phi.bc().type(b) == BndType::wall() ||
        phi.bc().type(b) == BndType::outlet() ) {

      Dir d = phi.bc().direction(b);

      if( d == Dir::imin() ) 
        for_vjk(phi.bc().at(b),j,k) 
         {fbnd[si()][j][k] += phi.bc().value(b) * dSx(Sign::neg(),si(),j,k);
          A.c[si()][j][k] -= A.w[si()][j][k];
          A.w[si()][j][k]  = 0.0;}

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {fbnd[ei()][j][k] += phi.bc().value(b) * dSx(Sign::pos(),ei(),j,k);
          A.c[ei()][j][k] -= A.e[ei()][j][k];
          A.e[ei()][j][k]  = 0.0;}

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {fbnd[i][sj()][k] += phi.bc().value(b) * dSy(Sign::neg(),i,sj(),k);
          A.c[i][sj()][k] -= A.s[i][sj()][k];
          A.s[i][sj()][k]  = 0.0;}

      if( d == Dir::jmax() ) 
        for_vik(phi.bc().at(b),i,k)
         {fbnd[i][ej()][k] += phi.bc().value(b) * dSy(Sign::pos(),i,ej(),k);
          A.c[i][ej()][k] -= A.n[i][ej()][k];
          A.n[i][ej()][k]  = 0.0;}

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {fbnd[i][j][sk()] += phi.bc().value(b) * dSz(Sign::neg(),i,j,sk());
          A.c[i][j][sk()] -= A.b[i][j][sk()];
          A.b[i][j][sk()]  = 0.0;}

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {fbnd[i][j][ek()] += phi.bc().value(b) * dSz(Sign::pos(),i,j,ek());
          A.c[i][j][ek()] -= A.t[i][j][ek()];
          A.t[i][j][ek()]  = 0.0;}

      if( d == Dir::ibody() ) {
        real area = 0.0;
        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

          if( phi.bc().at(b).contains_ijk(i,j,k) ) {

            area += dom->ibody().dS(cc);

            /* some cells are inside the immersed body 
               (there seems to be  too many checks, each
               face is checked twice, but it is needed) */
            if( dom->ibody().off(i,j,k) || dom->ibody().off(i-1,j,k) ) {
              A.c[i]  [j][k] -= A.w[i]  [j][k]; A.w[i]  [j][k] = 0.0;
              A.c[i-1][j][k] -= A.e[i-1][j][k]; A.e[i-1][j][k] = 0.0;
            }
            if( dom->ibody().off(i,j,k) || dom->ibody().off(i+1,j,k) ) {
              A.c[i]  [j][k] -= A.e[i]  [j][k]; A.e[i]  [j][k] = 0.0;
              A.c[i+1][j][k] -= A.w[i+1][j][k]; A.w[i+1][j][k] = 0.0;
            }
  
            if( dom->ibody().off(i,j,k) || dom->ibody().off(i,j-1,k) ) {
              A.c[i][j]  [k] -= A.s[i][j]  [k]; A.s[i][j]  [k] = 0.0;
              A.c[i][j-1][k] -= A.n[i][j-1][k]; A.n[i][j-1][k] = 0.0;
            }
            if( dom->ibody().off(i,j,k) || dom->ibody().off(i,j+1,k) ) {
              A.c[i][j]  [k] -= A.n[i][j]  [k]; A.n[i][j]  [k] = 0.0;
              A.c[i][j+1][k] -= A.s[i][j+1][k]; A.s[i][j+1][k] = 0.0;
            }
  
            if( dom->ibody().off(i,j,k) || dom->ibody().off(i,j,k-1) ) {
              A.c[i][j][k]   -= A.b[i][j][k];   A.b[i][j][k]   = 0.0;
              A.c[i][j][k-1] -= A.t[i][j][k-1]; A.t[i][j][k-1] = 0.0;
            }
            if( dom->ibody().off(i,j,k) || dom->ibody().off(i,j,k+1) ) {
              A.c[i][j][k]   -= A.t[i][j][k];   A.t[i][j][k]   = 0.0;
              A.c[i][j][k+1] -= A.b[i][j][k+1]; A.b[i][j][k+1] = 0.0;
            }

            /* compute source terms */
            if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) ) 
              fbnd[i][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            if( dom->ibody().on(i-1,j,k) && dom->ibody().off(i,j,k) ) 
              fbnd[i-1][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            
            if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) ) 
              fbnd[i][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            if( dom->ibody().on(i+1,j,k) && dom->ibody().off(i,j,k) ) 
              fbnd[i+1][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
    
            if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) ) 
              fbnd[i][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            if( dom->ibody().on(i,j-1,k) && dom->ibody().off(i,j,k) ) 
              fbnd[i][j-1][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            
            if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) ) 
              fbnd[i][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            if( dom->ibody().on(i,j+1,k) && dom->ibody().off(i,j,k) ) 
              fbnd[i][j+1][k] += phi.bc().value(b) * dom->ibody().dS(cc);
  
            if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) ) 
              fbnd[i][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            if( dom->ibody().on(i,j,k-1) && dom->ibody().off(i,j,k) ) 
              fbnd[i][j][k-1] += phi.bc().value(b) * dom->ibody().dS(cc);
          
            if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) ) 
              fbnd[i][j][k] += phi.bc().value(b) * dom->ibody().dS(cc);
            if( dom->ibody().on(i,j,k+1) && dom->ibody().off(i,j,k) ) 
              fbnd[i][j][k+1] += phi.bc().value(b) * dom->ibody().dS(cc);
          }
        }
      }
    }
  }

  /*----------------------+ 
  |  dirichlet and inlet  |
  +----------------------*/
  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type(b) == BndType::dirichlet() ||
        phi.bc().type(b) == BndType::inlet()     ||
        phi.bc().type(b) == BndType::insert() ) {

      Dir d = phi.bc().direction(b);

      if( d == Dir::imin() ) 
        for_vjk(phi.bc().at(b),j,k) 
         {const real value = phi[si()-1][j][k];
          fbnd[si()][j][k] += value * A.w[si()][j][k];
          A.w[si()][j][k]  = 0.0;}

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {const real value = phi[ei()+1][j][k];
          fbnd[ei()][j][k] += value * A.e[ei()][j][k];
          A.e[ei()][j][k]  = 0.0;}

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {const real value = phi[i][sj()-1][k];
          fbnd[i][sj()][k] += value * A.s[i][sj()][k];
          A.s[i][sj()][k]  = 0.0;}

      if( d == Dir::jmax() )  
        for_vik(phi.bc().at(b),i,k)
         {const real value = phi[i][ej()+1][k]; 
          fbnd[i][ej()][k] += value * A.n[i][ej()][k];
          A.n[i][ej()][k]  = 0.0;}

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {const real value = phi[i][j][sk()-1];
          fbnd[i][j][sk()] += value * A.b[i][j][sk()];
          A.b[i][j][sk()]  = 0.0;}

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {const real value = phi[i][j][ek()+1];
          fbnd[i][j][ek()] += value * A.t[i][j][ek()];
          A.t[i][j][ek()]  = 0.0;}

      if( d == Dir::ibody() ) {

        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) ) 
           {fbnd[i][j][k] += phi.bc().value(b) * A.w[i][j][k];
            A.w [i][j][k] = 0.0;}
          if( dom->ibody().on(i-1,j,k) && dom->ibody().off(i,j,k) ) 
           {fbnd[i-1][j][k] += phi.bc().value(b) * A.e[i-1][j][k];
            A.e [i-1][j][k] = 0.0;}
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) ) 
           {fbnd[i][j][k] += phi.bc().value(b) * A.e[i][j][k];
            A.e [i][j][k] = 0.0;}
          if( dom->ibody().on(i+1,j,k) && dom->ibody().off(i,j,k) ) 
           {fbnd[i+1][j][k] += phi.bc().value(b) * A.w[i+1][j][k];
            A.w [i+1][j][k] = 0.0;}
          
  
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) ) 
           {fbnd[i][j][k] += phi.bc().value(b) * A.s[i][j][k];
            A.s [i][j][k] = 0.0;}
          if( dom->ibody().on(i,j-1,k) && dom->ibody().off(i,j,k) ) 
           {fbnd[i][j-1][k] += phi.bc().value(b) * A.n[i][j-1][k];
            A.n [i][j-1][k] = 0.0;}
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) ) 
           {fbnd[i][j][k] += phi.bc().value(b) * A.n[i][j][k];
            A.n [i][j][k] = 0.0;}
          if( dom->ibody().on(i,j+1,k) && dom->ibody().off(i,j,k) ) 
           {fbnd[i][j+1][k] += phi.bc().value(b) * A.s[i][j+1][k];
            A.s [i][j+1][k] = 0.0;}
          
  
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) ) 
           {fbnd[i][j][k] += phi.bc().value(b) * A.b[i][j][k];
            A.b [i][j][k] = 0.0;}
          if( dom->ibody().on(i,j,k-1) && dom->ibody().off(i,j,k) ) 
           {fbnd[i][j][k-1] += phi.bc().value(b) * A.t[i][j][k-1];
            A.t [i][j][k-1] = 0.0;}
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) ) 
           {fbnd[i][j][k] += phi.bc().value(b) * A.t[i][j][k];
            A.t [i][j][k] = 0.0;}
          if( dom->ibody().on(i,j,k+1) && dom->ibody().off(i,j,k) ) 
           {fbnd[i][j][k+1] += phi.bc().value(b) * A.b[i][j][k+1];
            A.b [i][j][k+1] = 0.0;}
          
        }
      }
    }
  }

  /*-------------+ 
  |  convective  |
  +-------------*/
  /* special treatment is needed here, since the boundary values must be 
     updated using heat transfer coefficient (coeff) and infinite
     temperatre. values at the boundary must be updated here too, since 
     these data (coeff and infinite) are not available in the Scalar class 
     itsel. */
  real tsc = diff_ts.N();

  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type(b) == BndType::convective() ) {

      Dir d = phi.bc().direction(b);

      if( d == Dir::imin() ) 
        for_vjk(phi.bc().at(b),j,k) 
         {const real value  = phi.bc().value(b, Comp::inf());
          const real coeff  = phi.bc().value(b, Comp::coefficient());
          const real lambda = f_prop -> value(si(),j,k);
          const real area   = dSx(Sign::neg(),si(),j,k);
          const real delta  = dxw(si());
          const real factor = (lambda/delta)/(coeff+lambda/delta)*coeff*area;
          fbnd[si()][j][k] += value * factor;
          A.c[si()][j][k]  += factor / tsc;
          phi[si()-1][j][k] = (value*coeff + phi[si()][j][k]*lambda/delta)
                            / (coeff + lambda/delta); 
          A.w[si()][j][k]  = 0.0;}

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {const real value  = phi.bc().value(b, Comp::inf());
          const real coeff  = phi.bc().value(b, Comp::coefficient());
          const real lambda = f_prop -> value(si(),j,k);
          const real area   = dSx(Sign::pos(),ei(),j,k);
          const real delta  = dxe(ei());
          const real factor = (lambda/delta)/(coeff+lambda/delta)*coeff*area;
          fbnd[ei()][j][k] += value * factor;
          A.c[ei()][j][k]  += factor / tsc;
          phi[ei()+1][j][k] = (value*coeff + phi[ei()][j][k]*lambda/delta)
                            / (coeff + lambda/delta); 
          A.e[ei()][j][k]  = 0.0;}

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {const real value  = phi.bc().value(b, Comp::inf());
          const real coeff  = phi.bc().value(b, Comp::coefficient());
          const real lambda = f_prop -> value(i,sj(),k);
          const real area   = dSy(Sign::neg(),i,sj(),k);
          const real delta  = dys(sj());
          const real factor = (lambda/delta)/(coeff+lambda/delta)*coeff*area;
          fbnd[i][sj()][k] += value * factor;
          A.c[i][sj()][k]  += factor / tsc;
          phi[i][sj()-1][k] = (value*coeff + phi[i][sj()][k]*lambda/delta)
                            / (coeff + lambda/delta); 
          A.s[i][sj()][k]  = 0.0;}

      if( d == Dir::jmax() )  
        for_vik(phi.bc().at(b),i,k)
         {const real value  = phi.bc().value(b, Comp::inf());
          const real coeff  = phi.bc().value(b, Comp::coefficient());
          const real lambda = f_prop -> value(i,ej(),k);
          const real area   = dSy(Sign::pos(),i,ej(),k);
          const real delta  = dyn(ej());
          const real factor = (lambda/delta)/(coeff+lambda/delta)*coeff*area;
          fbnd[i][ej()][k] += value * factor;
          A.c[i][ej()][k]  += factor / tsc;
          phi[i][ej()+1][k] = (value*coeff + phi[i][ej()][k]*lambda/delta)
                            / (coeff + lambda/delta); 
          A.n[i][ej()][k]  = 0.0;}

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {const real value  = phi.bc().value(b, Comp::inf());
          const real coeff  = phi.bc().value(b, Comp::coefficient());
          const real lambda = f_prop -> value(i,j,sk());
          const real area   = dSz(Sign::neg(),i,j,sk());
          const real delta  = dzb(sk());
          const real factor = (lambda/delta)/(coeff+lambda/delta)*coeff*area;
          fbnd[i][j][sk()] += value * factor;
          A.c[i][j][sk()]  += factor / tsc;
          phi[i][j][sk()-1] = (value*coeff + phi[i][j][sk()]*lambda/delta)
                            / (coeff + lambda/delta); 
          A.b[i][j][sk()]  = 0.0;}

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {const real value  = phi.bc().value(b, Comp::inf());
          const real coeff  = phi.bc().value(b, Comp::coefficient());
          const real lambda = f_prop -> value(i,j,ek());
          const real area   = dSz(Sign::pos(),i,j,ek());
          const real delta  = dzt(ek());
          const real factor = (lambda/delta)/(coeff+lambda/delta)*coeff*area;
          fbnd[i][j][ek()] += value * factor;
          A.c[i][j][ek()]  += factor / tsc;
          phi[i][j][ek()+1] = (value*coeff + phi[i][j][ek()]*lambda/delta)
                            / (coeff + lambda/delta); 
          A.t[i][j][ek()]  = 0.0;}

      /* oh well, what to do with immersed body? */
      if( d == Dir::ibody() ) {

      }
    }
  }

  A.c.exchange();
  A.w.exchange();
  A.e.exchange();
  A.s.exchange();
  A.n.exchange();
  A.b.exchange();
  A.t.exchange();
  A.ci.exchange();

}
