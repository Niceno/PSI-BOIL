#include "momentum.h"

/******************************************************************************/
real Momentum::volf_bct(const BndType & bc_type,
                        real * Ax, real * Ay, real * Az) const {
/*----------------------------------------------------------+ 
|  computes volume flux [m^3/s] at specified b.c. type      |
+- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
|  no density is used - the flow is incompressible          |
+----------------------------------------------------------*/

  int  m;
  real a_x, a_y, a_z;
  real blk = 0.0;
  real arex = 0.0;
  real arey = 0.0;
  real arez = 0.0;

  /*-------------------------+ 
  |  through all boundaries  |
  +-------------------------*/
  for_m(m)
  for( int b=0; b<u.bc(m).count(); b++ ) {

    if( u.bc(m).type(b) == bc_type ) {

      Dir d = u.bc(m).direction(b);

      /*------------------------------+
      |  through velocity components  |
      +------------------------------*/

      /*--------------+
      |  imin & imax  | 
      +--------------*/
      if(m==Comp::u()) {

        int imin = (d==Dir::imin()); 
        int imax = (d==Dir::imax()); 
        real volf_x = 0.0;
        real A_x    = 0.0;

        if( imin || imax ) {

          if(imin && dom->coord(Comp::i())==0) {
            const int i=si(m);
            for_vjk(u.bc(m).at(b),j,k) {
              if( dom->ibody().on(m, i, j, k) ) {
                a_x = dSx(m,i,j,k); 
                a_x *= dom->ibody().fSw(m,i,j,k);              
                volf_x += a_x * u[m][i-1][j][k];
                A_x    += a_x;
              }
            }
            //std::cout<<"volf_bct imin"<<A_x<<"\n";
          }
          if(imax && dom->coord(Comp::i())==dom->dim(Comp::i())-1) {
            const int i=ei(m);
            for_vjk(u.bc(m).at(b),j,k) {
              if( dom->ibody().on(m, i, j, k) ) {
                a_x = dSx(m,i,j,k); 
                a_x *= dom->ibody().fSe(m,i,j,k);              
                volf_x -= a_x * u[m][i+1][j][k];
                A_x    += a_x;
                //boil::oout<<i<<" "<<j<<" "<<k<<" "<<a_x<<boil::endl;
              }
            }
            //std::cout<<"volf_bct imax"<<A_x<<"\n";
            //exit(0);
          }

        } /* if(imin || imax) */

        boil::cart.sum_int(&imin);
        boil::cart.sum_int(&imax);

        if( imin || imax ) {
          boil::cart.sum_real(&volf_x);
          boil::cart.sum_real(&A_x);
          blk += volf_x; 
          arex += A_x;   
          if(Ax != NULL) *Ax = arex;
        } /* if(imin || imax) */

      } /* i block */

      /*--------------+
      |  jmin & jmax  | 
      +--------------*/
      if(m==Comp::v()) {
        
        int  jmin = (d==Dir::jmin()); 
        int  jmax = (d==Dir::jmax()); 
        real volf_y = 0.0;
        real A_y    = 0.0;
  
        if( jmin || jmax ) {

          if(jmin && dom->coord(Comp::j())==0) {
            const int j=sj(m);
            for_vik(u.bc(m).at(b),i,k) {
              if( dom->ibody().on(m, i, j, k) ) {
                a_y = dSy(m,i,j,k); 
                a_y *= dom->ibody().fSs(m,i,j,k);
                volf_y += a_y * u[m][i][j-1][k];
                A_y    += a_y;
              }
            }
            //std::cout<<"volf_bct jmin"<<A_y<<"\n";
          }
          if(jmax && dom->coord(Comp::j())==dom->dim(Comp::j())-1) {
            const int j=ej(m);
            for_vik(u.bc(m).at(b),i,k) {
              if( dom->ibody().on(m, i, j, k) ) {
                a_y = dSy(m,i,j,k); 
                a_y *= dom->ibody().fSn(m,i,j,k);
                volf_y -= a_y * u[m][i][j+1][k];
                A_y    += a_y;
              }
            }
            //std::cout<<"volf_bct jmax"<<A_y<<"\n";
          }

        } /* if(jmin || jmax) */

        boil::cart.sum_int(&jmin);
        boil::cart.sum_int(&jmax);

        if( jmin || jmax ) {
          boil::cart.sum_real(&volf_y);
          boil::cart.sum_real(&A_y);
          blk += volf_y;
          arey += A_y;
          if(Ay != NULL) *Ay = arey;
        } /* if(jmin || jmax) */

      } /* j block */

      /*--------------+
      |  kmin & kmax  | 
      +--------------*/
      if(m==Comp::w()) {

        int  kmin = (d==Dir::kmin()); 
        int  kmax = (d==Dir::kmax()); 
        real volf_z = 0.0;
        real A_z    = 0.0;
  
        if( kmin || kmax ) {

          if(kmin && dom->coord(Comp::k())==0) {
            const int k=sk(m);
            for_vij(u.bc(m).at(b),i,j) {
              if( dom->ibody().on(m, i, j, k) ) {
                a_z = dSz(m,i,j,k); 
                a_z *= dom->ibody().fSb(m,i,j,k);              
                volf_z += a_z * u[m][i][j][k-1];
                //volf_z += a_z * u[m][i][j][k+1];
                A_z    += a_z;
              }
            }
          }
          if(kmax && dom->coord(Comp::k())==dom->dim(Comp::k())-1) {
            const int k=ek(m);
            for_vij(u.bc(m).at(b),i,j) {
              if( dom->ibody().on(m, i, j, k) ) {
                a_z = dSz(m,i,j,k); 
                a_z *= dom->ibody().fSt(m,i,j,k);              
                volf_z -= a_z * u[m][i][j][k+1];
                //volf_z -= a_z * u[m][i][j][k-1];
                A_z    += a_z;
              }
            }
          }
        } /* if(kmin || kmax) */

        boil::cart.sum_int(&kmin);
        boil::cart.sum_int(&kmax);

        if( kmin || kmax ) {
          boil::cart.sum_real(&volf_z);
          boil::cart.sum_real(&A_z);
          blk += volf_z; 
          arez += A_z;   
          if(Az != NULL) *Az = arez;
        } /* if(kmin || kmax) */

      } /* k block */

    } /* if desired type */

  } /* b */

  //assert(are > 0.0);

  return blk;
}
