#include "vector.h"

/******************************************************************************/
real Vector::bnd_flow(const BndType & bc_type,
                      real * Ax, real * Ay, real * Az) const {
/*----------------------------------------------------------+ 
|  computes flow [units * m^2] at specified b.c. type       |
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
  for( int b=0; b<bc(m).count(); b++ ) {

    if( bc(m).type(b) == bc_type ) {

      Dir d = bc(m).direction(b);

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

          if(imin && domain()->coord(Comp::i())==0) {
            const int i=si(m);
            for_vjk(bc(m).at(b),j,k) {
              if( domain()->ibody().on(m, i, j, k) ) {
                //a_x = dSx(m,i,j,k); 
                a_x = domain()->dSx(Sign::neg(),i,j,k); 
                a_x *= domain()->ibody().fSw(m,i,j,k);              
                volf_x += a_x * vec[m][i-1][j][k];
                A_x    += a_x;
              }
            }
            //std::cout<<"volf_bct imin"<<A_x<<"\n";
          }
          if(imax && domain()->coord(Comp::i())==domain()->dim(Comp::i())-1) {
            const int i=ei(m);
            for_vjk(bc(m).at(b),j,k) {
              if( domain()->ibody().on(m, i, j, k) ) {
                //a_x = dSx(m,i,j,k); 
                a_x = domain()->dSx(Sign::pos(),i,j,k); 
                a_x *= domain()->ibody().fSe(m,i,j,k);              
                volf_x -= a_x * vec[m][i+1][j][k];
                A_x    += a_x;
                //boil::oout<<i<<" "<<j<<" "<<k<<" "<<a_x<<" "<<vec[m][i+1][j][k]<<boil::endl;
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

          if(jmin && domain()->coord(Comp::j())==0) {
            const int j=sj(m);
            for_vik(bc(m).at(b),i,k) {
              if( domain()->ibody().on(m, i, j, k) ) {
                //a_y = dSy(m,i,j,k); 
                a_y = domain()->dSy(Sign::neg(),i,j,k);
                a_y *= domain()->ibody().fSs(m,i,j,k);
                volf_y += a_y * vec[m][i][j-1][k];
                A_y    += a_y;
              }
            }
            //std::cout<<"volf_bct jmin"<<A_y<<"\n";
          }
          if(jmax && domain()->coord(Comp::j())==domain()->dim(Comp::j())-1) {
            const int j=ej(m);
            for_vik(bc(m).at(b),i,k) {
              if( domain()->ibody().on(m, i, j, k) ) {
                //a_y = dSy(m,i,j,k); 
                a_y = domain()->dSy(Sign::pos(),i,j,k);
                a_y *= domain()->ibody().fSn(m,i,j,k);
                volf_y -= a_y * vec[m][i][j+1][k];
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

          if(kmin && domain()->coord(Comp::k())==0) {
            const int k=sk(m);
            for_vij(bc(m).at(b),i,j) {
              if( domain()->ibody().on(m, i, j, k) ) {
                //a_z = dSz(m,i,j,k); 
                a_z = domain()->dSz(Sign::neg(),i,j,k);
                a_z *= domain()->ibody().fSb(m,i,j,k);              
                volf_z += a_z * vec[m][i][j][k-1];
                //volf_z += a_z * vec[m][i][j][k+1];
                A_z    += a_z;
              }
            }
          }
          if(kmax && domain()->coord(Comp::k())==domain()->dim(Comp::k())-1) {
            const int k=ek(m);
            for_vij(bc(m).at(b),i,j) {
              if( domain()->ibody().on(m, i, j, k) ) {
                //a_z = dSz(m,i,j,k); 
                a_z = domain()->dSz(Sign::pos(),i,j,k);
                a_z *= domain()->ibody().fSt(m,i,j,k);              
                volf_z -= a_z * vec[m][i][j][k+1];
                //volf_z -= a_z * vec[m][i][j][k-1];
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
