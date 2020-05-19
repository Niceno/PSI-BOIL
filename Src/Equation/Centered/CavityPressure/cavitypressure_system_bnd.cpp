#include "cavitypressure.h"

/***************************************************************************//**
*  \brief Corrects system matrix \f$ [A] \f$ at boundaries.
*******************************************************************************/
void CavityPressure::create_system_bnd() {

  /*--------------------------------------+
  |  initialize, that is quite important  |
  +--------------------------------------*/
  fbnd = 0.0;
  phi.bnd_update();

  std::function<bool(const int,const int,const int)> valid_cell
    = [this](const int i, const int j, const int k)
    { return dom->ibody().on_p(i,j,k)
             && !in_gas(i,j,k); };

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
         {//fbnd[si()][j][k] += phi.bc().value(b) * dSx(Sign::neg(),si(),j,k);
          if(valid_cell(si(),j,k)&&!interface(Sign::neg(),Comp::i(),si(),j,k)) {
            A.c[si()][j][k] -= A.w[si()][j][k];
            A.w[si()][j][k]  = 0.0;
          }
         }

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {//fbnd[ei()][j][k] += phi.bc().value(b) * dSx(Sign::pos(),ei(),j,k);
          if(valid_cell(ei(),j,k)&&!interface(Sign::pos(),Comp::i(),ei(),j,k)) {
            A.c[ei()][j][k] -= A.e[ei()][j][k];
            A.e[ei()][j][k]  = 0.0;
          }
         }

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {//fbnd[i][sj()][k] += phi.bc().value(b) * dSy(Sign::neg(),i,sj(),k);
          if(valid_cell(i,sj(),k)&&!interface(Sign::neg(),Comp::j(),i,sj(),k)) {
            A.c[i][sj()][k] -= A.s[i][sj()][k];
            A.s[i][sj()][k]  = 0.0;
          }
         }

      if( d == Dir::jmax() ) 
        for_vik(phi.bc().at(b),i,k)
         {//fbnd[i][ej()][k] += phi.bc().value(b) * dSy(Sign::pos(),i,ej(),k);
          if(valid_cell(i,ej(),k)&&!interface(Sign::pos(),Comp::j(),i,ej(),k)) {
            A.c[i][ej()][k] -= A.n[i][ej()][k];
            A.n[i][ej()][k]  = 0.0;
          }
         }

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {//fbnd[i][j][sk()] += phi.bc().value(b) * dSz(Sign::neg(),i,j,sk());
          if(valid_cell(i,j,sk())&&!interface(Sign::neg(),Comp::k(),i,j,sk())) {
            A.c[i][j][sk()] -= A.b[i][j][sk()];
            A.b[i][j][sk()]  = 0.0;
          }
         }

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {//fbnd[i][j][ek()] += phi.bc().value(b) * dSz(Sign::pos(),i,j,ek());
          if(valid_cell(i,j,ek())&&!interface(Sign::pos(),Comp::k(),i,j,ek())) {
            A.c[i][j][ek()] -= A.t[i][j][ek()];
            A.t[i][j][ek()]  = 0.0;
          }
         }

      /* this part should be deprecated! */
      if( d == Dir::ibody() ) {
        boil::oout << "CavityPressure::system_bnd: Underdevelopment!"<<boil::endl;
        exit(0);
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
          if(valid_cell(si(),j,k)&&!interface(Sign::neg(),Comp::i(),si(),j,k)) {
            const real value = phi[si()-1][j][k];
            fbnd[si()][j][k] += value * A.w[si()][j][k];
            A.w[si()][j][k]  = 0.0;
          }
         }

      if( d == Dir::imax() ) 
        for_vjk(phi.bc().at(b),j,k)
         {
          if(valid_cell(ei(),j,k)&&!interface(Sign::pos(),Comp::i(),ei(),j,k)) {
            const real value = phi[ei()+1][j][k];
            fbnd[ei()][j][k] += value * A.e[ei()][j][k];
            boil::oout<<k<<" "<<value<<" "<<fbnd[ei()][j][k]<<" "<<A.e[ei()][j][k]<<boil::endl;
            A.e[ei()][j][k]  = 0.0;
          }
         }

      if( d == Dir::jmin() ) 
        for_vik(phi.bc().at(b),i,k)
         {
          if(valid_cell(i,sj(),k)&&!interface(Sign::neg(),Comp::j(),i,sj(),k)) {
            const real value = phi[i][sj()-1][k];
            fbnd[i][sj()][k] += value * A.s[i][sj()][k];
            A.s[i][sj()][k]  = 0.0;
          }
         }

      if( d == Dir::jmax() )  
        for_vik(phi.bc().at(b),i,k)
         {
          if(valid_cell(i,ej(),k)&&!interface(Sign::pos(),Comp::j(),i,ej(),k)) {
            const real value = phi[i][ej()+1][k];
            fbnd[i][ej()][k] += value * A.n[i][ej()][k];
            A.n[i][ej()][k]  = 0.0;
          }
         }

      if( d == Dir::kmin() ) 
        for_vij(phi.bc().at(b),i,j)
         {
          if(valid_cell(i,j,sk())&&!interface(Sign::neg(),Comp::k(),i,j,sk())) {
            const real value = phi[i][j][sk()-1];
            fbnd[i][j][sk()] += value * A.b[i][j][sk()];
            A.b[i][j][sk()]  = 0.0;
          }
         }

      if( d == Dir::kmax() ) 
        for_vij(phi.bc().at(b),i,j)
         {
          if(valid_cell(i,j,ek())&&!interface(Sign::pos(),Comp::k(),i,j,ek())) {
            const real value = phi[i][j][ek()+1];
            fbnd[i][j][ek()] += value * A.t[i][j][ek()];
            A.t[i][j][ek()]  = 0.0;
          }
         }

      /* this part should be deprecated! */
      if( d == Dir::ibody() ) {
        boil::oout << "CavityPressure::system_bnd: Underdevelopment!"<<boil::endl;
        exit(0);
      } /* dir = ibody */
    }
  }

  /*----------------+ 
  |  immersed body  |
  +----------------*/
  if(dom->ibody().nccells() > 0) {
    for(int cc=0; cc<dom->ibody().nccells(); cc++) {
      int i,j,k;
      /* i,j,k is cell adjacent to solid */
      dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

      if(valid_cell(i,j,k)) {
        if(dom->ibody().off_p(i-1,j,k)&&!interface(Sign::neg(),Comp::i(),i,j,k))
          {A.c[i][j][k] -= A.w[i][j][k]; A.w[i][j][k] = 0.0;}
        if(dom->ibody().off_p(i+1,j,k)&&!interface(Sign::pos(),Comp::i(),i,j,k))
          {A.c[i][j][k] -= A.e[i][j][k]; A.e[i][j][k] = 0.0;}
        if(dom->ibody().off_p(i,j-1,k)&&!interface(Sign::neg(),Comp::j(),i,j,k))
          {A.s[i][j][k] -= A.s[i][j][k]; A.s[i][j][k] = 0.0;}
        if(dom->ibody().off_p(i,j+1,k)&&!interface(Sign::pos(),Comp::j(),i,j,k))
          {A.n[i][j][k] -= A.n[i][j][k]; A.n[i][j][k] = 0.0;}
        if(dom->ibody().off_p(i,j,k-1)&&!interface(Sign::neg(),Comp::k(),i,j,k))
          {A.b[i][j][k] -= A.b[i][j][k]; A.b[i][j][k] = 0.0;}
        if(dom->ibody().off_p(i,j,k+1)&&!interface(Sign::pos(),Comp::k(),i,j,k))
          {A.t[i][j][k] -= A.t[i][j][k]; A.t[i][j][k] = 0.0;}
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

}	

