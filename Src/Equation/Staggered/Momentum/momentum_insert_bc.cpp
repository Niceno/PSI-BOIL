#include "momentum.h"

/******************************************************************************/
void Momentum::insert_bc() {
/*----------------------------+ 
|  insert boundary condition  |
+----------------------------*/

  Formula F;

  /* check boundary conditions */
  for_m(m)
    for( int b=0; b<u.bc(m).count(); b++ ) {
 
      if( u.bc(m).type(b) == BndType::dirichlet() ||
          u.bc(m).type(b) == BndType::neumann() ) {
        boil::oout << "Unsupported b.c. for momentum at " << b << boil::endl; 
        exit(0);
      }
    }

  /* set-up the symmetry in a better way */
  for_m(m)
    for( int b=0; b<u.bc(m).count(); b++ ) {
      if( u.bc(m).type(b) == BndType::symmetry() ) {
        Dir d = u.bc(m).direction(b);
        if(
           m == Comp::u() && ( d == Dir::imin() || d == Dir::imax()) ||
           m == Comp::v() && ( d == Dir::jmin() || d == Dir::jmax()) ||
           m == Comp::w() && ( d == Dir::kmin() || d == Dir::kmax()) 
          ) {
              u.bc(m).type(b) = BndType::wall();
              boil::oout << "Adjusting symmetry b.c. for momentum at " << b
                         << boil::endl;
        } else {  // 28.11.2013 Yohei
          if( !u.bc(m).type_decomp(b) ) {
          for_vijk(u.bc(m).at(b), i,j,k) {
            int iof=0, jof=0, kof=0;
            if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
            if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
            if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;
            u[m][i][j][k] = u[m][i+iof][j+jof][k+kof];
          }
          }
        }
      } /* symmetry */
    } /* b */

  /* inlet and wall */
  for_m(m)
    for( int b=0; b<u.bc(m).count(); b++ ) {

      if( u.bc(m).exists(b) ) {

        if( u.bc(m).type(b) == BndType::inlet() ||
            u.bc(m).type(b) == BndType::wall() ) {

          /* formula is defined */
          if( u.bc(m).formula(b,m) ) {

            for_vijk(u.bc(m).at(b), i,j,k) {
              std::stringstream x, y, z, f;
              x << "x=" << u.xc(m,i); F.evaluate(x);
              y << "y=" << u.yc(m,j); F.evaluate(y);
              z << "z=" << u.zc(m,k); F.evaluate(z);
              f << u.bc(m).formula(b,m);
              u[m][i][j][k] = F.evaluate(f);
            }

          /* formula is not defined */
          } else {
            for_vijk(u.bc(m).at(b), i,j,k){
              u[m][i][j][k] = u.bc(m).value(b,m);
            }
          }

        } /* if wall or inlet */
      } /* if boundary exists */
    } /* m & b */
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_insert_bc.cpp,v 1.13 2014/02/04 13:36:21 sato Exp $'/
+-----------------------------------------------------------------------------*/
