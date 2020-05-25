#include "vector.h"

/******************************************************************************/
void Vector::bnd_update_nooutlet() {
/*----------------------------+ 
|  insert boundary condition  |
+----------------------------*/

  Formula F;

  /* check boundary conditions */
  for_m(m)
    for( int b=0; b<bc(m).count(); b++ ) {
 
      if(   bc(m).type(b) == BndType::dirichlet()
        // || bc(m).type(b) == BndType::neumann() 
        ) {
        boil::oout << "Unsupported b.c. for momentum at " << b << boil::endl; 
        exit(0);
      }
    }

  /* set-up the symmetry in a better way */
  for_m(m)
    for( int b=0; b<bc(m).count(); b++ ) {
      if(   bc(m).type(b) == BndType::symmetry()
         || bc(m).type(b) == BndType::neumann() ) {
        Dir d = bc(m).direction(b);
        if(
           m == Comp::u() && ( d == Dir::imin() || d == Dir::imax()) ||
           m == Comp::v() && ( d == Dir::jmin() || d == Dir::jmax()) ||
           m == Comp::w() && ( d == Dir::kmin() || d == Dir::kmax()) 
          ) {
              bc(m).type(b) = BndType::wall();
              boil::oout << "Adjusting symmetry b.c. for momentum at " << b
                         << boil::endl;
        } else {  // 28.11.2013 Yohei
          if( !bc(m).type_decomp(b) ) {
          for_vijk(bc(m).at(b), i,j,k) {
            int iof=0, jof=0, kof=0;
            if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
            if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
            if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;
            vec[m][i][j][k] = vec[m][i+iof][j+jof][k+kof];
          }
          }
        }
      } /* symmetry */
    } /* b */

  /* inlet and wall */
  for_m(m)
    for( int b=0; b<bc(m).count(); b++ ) {

      if( bc(m).exists(b) ) {

        if( bc(m).type(b) == BndType::inlet() ||
            bc(m).type(b) == BndType::wall() ) {

          /* formula is defined */
          if( bc(m).formula(b,m) ) {

            for_vijk(bc(m).at(b), i,j,k) {
              std::stringstream x, y, z, f;
              x << "x=" << xc(m,i); F.evaluate(x);
              y << "y=" << yc(m,j); F.evaluate(y);
              z << "z=" << zc(m,k); F.evaluate(z);
              f << bc(m).formula(b,m);
              vec[m][i][j][k] = F.evaluate(f);
            }

          /* formula is not defined */
          } else {
            for_vijk(bc(m).at(b), i,j,k){
              vec[m][i][j][k] = bc(m).value(b,m);
            }
          }

        } /* if wall or inlet */
      } /* if boundary exists */
    } /* m & b */

#if 1
  /* pseudo b.c. */
  for_m(m)
    for( int b=0; b<bc(m).count(); b++ ) {
      if( bc(m).type(b) == BndType::pseudo() ) {
        Dir d = bc(m).direction(b);
        if(
           m == Comp::u() && ( d != Dir::imin() && d != Dir::imax()) ||
           m == Comp::v() && ( d != Dir::jmin() && d != Dir::jmax()) ||
           m == Comp::w() && ( d != Dir::kmin() && d != Dir::kmax())
          ) {
          int ii(0), jj(0), kk(0);
          if(d == Dir::imin()) ii =  1;
          if(d == Dir::imax()) ii = -1;
          if(d == Dir::jmin()) jj =  1;
          if(d == Dir::jmax()) jj = -1;
          if(d == Dir::kmin()) kk =  1;
          if(d == Dir::kmax()) kk = -1;
          for_vijk(bc(m).at(b), i,j,k) {
            vec[m][i][j][k] = vec[m][i+ii][j+jj][k+kk];
          }
        }  
      }
    }
#endif

  return;
  
}
