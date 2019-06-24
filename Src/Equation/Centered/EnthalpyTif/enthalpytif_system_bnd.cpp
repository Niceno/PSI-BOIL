#include "enthalpytif.h"

/***************************************************************************//**
*  \brief Corrects system matrix \f$ [A] \f$ at boundaries.
*******************************************************************************/
void EnthalpyTIF::create_system_bnd() {

  /*--------------------------------------+
  |  initialize, that is quite important  |
  +--------------------------------------*/
  fbnd = 0.0;

  /*----------------------+ 
  |  symmetry and outlet  |
  +----------------------*/
  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type(b) == BndType::neumann() ||
        phi.bc().type(b) == BndType::symmetry()  ||
        phi.bc().type(b) == BndType::pseudo() ||
        phi.bc().type(b) == BndType::outlet() ) {

      Dir d = phi.bc().direction(b);

      if( d == Dir::imin() ) 
        for_vjk(phi.bc().at(b),j,k) 
         {//fbnd[si()][j][k] += phi.bc().value(b) * dSx(si(),j,k);
          if(!fs||!Interface(-1,Comp::i(),si(),j,k)) {
            A.c[si()][j][k] -= A.w[si()][j][k];
            A.w[si()][j][k]  = 0.0;
          }
         }

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {//fbnd[ei()][j][k] += phi.bc().value(b) * dSx(ei(),j,k);
          if(!fs||!Interface(+1,Comp::i(),ei(),j,k)) {
            A.c[ei()][j][k] -= A.e[ei()][j][k];
            A.e[ei()][j][k]  = 0.0;
          }
         }

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {//fbnd[i][sj()][k] += phi.bc().value(b) * dSy(i,sj(),k);
          if(!fs||!Interface(-1,Comp::j(),i,sj(),k)) {
            A.c[i][sj()][k] -= A.s[i][sj()][k];
            A.s[i][sj()][k]  = 0.0;
          }
         }

      if( d == Dir::jmax() ) 
        for_vik(phi.bc().at(b),i,k)
         {//fbnd[i][ej()][k] += phi.bc().value(b) * dSy(i,ej(),k);
          if(!fs||!Interface(+1,Comp::j(),i,ej(),k)) {
            A.c[i][ej()][k] -= A.n[i][ej()][k];
            A.n[i][ej()][k]  = 0.0;
          }
         }

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {//fbnd[i][j][sk()] += phi.bc().value(b) * dSz(i,j,sk());
          if(!fs||!Interface(-1,Comp::k(),i,j,sk())) {
            A.c[i][j][sk()] -= A.b[i][j][sk()];
            A.b[i][j][sk()]  = 0.0;
          }
         }

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {//fbnd[i][j][ek()] += phi.bc().value(b) * dSz(i,j,ek());
          if(!fs||!Interface(+1,Comp::k(),i,j,ek())) {
            A.c[i][j][ek()] -= A.t[i][j][ek()];
            A.t[i][j][ek()]  = 0.0;
          }
         }

      /* this part should be deprecated! */
      if( d == Dir::ibody() ) {
        boil::oout << "ETIF::system_bnd: Underdevelopment!"<<boil::endl;
        exit(0);
        #if 0
        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

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
        } /* loop over ibody cells */
        #endif
      } /* dir = ibody */
    } /* if condition */
  } /* loop over bc directions */

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
         {
          if(!fs||!Interface(-1,Comp::i(),si(),j,k)) {
            const real value = phi[si()-1][j][k];
            fbnd[si()][j][k] += value * A.w[si()][j][k];
            A.w[si()][j][k]  = 0.0;
          }
         }

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {
          if(!fs||!Interface(+1,Comp::i(),ei(),j,k)) {
            const real value = phi[ei()+1][j][k];
            fbnd[ei()][j][k] += value * A.e[ei()][j][k];
            A.e[ei()][j][k]  = 0.0;
          }
         }

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {
          if(!fs||!Interface(-1,Comp::j(),i,sj(),k)) {
            const real value = phi[i][sj()-1][k];
            fbnd[i][sj()][k] += value * A.s[i][sj()][k];
            A.s[i][sj()][k]  = 0.0;
          }
         }

      if( d == Dir::jmax() )  
        for_vik(phi.bc().at(b),i,k)
         {
          if(!fs||!Interface(+1,Comp::j(),i,ej(),k)) {
            const real value = phi[i][ej()+1][k];
            fbnd[i][ej()][k] += value * A.n[i][ej()][k];
            A.n[i][ej()][k]  = 0.0;
          }
         }

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {
          if(!fs||!Interface(-1,Comp::k(),i,j,sk())) {
            const real value = phi[i][j][sk()-1];
            fbnd[i][j][sk()] += value * A.b[i][j][sk()];
            A.b[i][j][sk()]  = 0.0;
          }
         }

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {
          if(!fs||!Interface(+1,Comp::k(),i,j,ek())) {
            const real value = phi[i][j][ek()+1];
            fbnd[i][j][ek()] += value * A.t[i][j][ek()];
            A.t[i][j][ek()]  = 0.0;
          }
         }

      /* this part should be deprecated! */
      if( d == Dir::ibody() ) {
        boil::oout << "ETIF::system_bnd: Underdevelopment!"<<boil::endl;
        boil::oout << "ibody without conduction in solid." <<boil::endl;
        exit(0);
        #if 0
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
        #endif
      } /* dir = ibody */
    }
  }

  A.c.exchange();
  A.w.exchange();
  A.e.exchange();
  A.s.exchange();
  A.n.exchange();
  A.b.exchange();
  A.t.exchange();

}	

