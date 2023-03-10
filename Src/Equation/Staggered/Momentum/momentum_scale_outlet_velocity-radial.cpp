#include "momentum.h"
//#define TWODXY
//#define TWODXZ
#define AXISYM
#define SYM

/******************************************************************************/
void Momentum::scale_outlet_velocity(const real ubo, const real ratio) {

#ifdef TWODXY
  #ifdef SYM
  const real prefact = 2.0/4.0*boil::pi;
  #else
  const real prefact = 2.0*boil::pi;
  #endif
#elif defined TWODXZ
  #ifdef SYM
  const real prefact = 2.0/4.0*boil::pi;
  #else
  const real prefact = 2.0*boil::pi;
  #endif
#elif defined AXISYM
  #ifdef SYM
  const real prefact = 1.0;
  #else
  const real prefact = 2.0;
  #endif
#else
  #ifdef SYM
  const real prefact = 4.0/8.0*boil::pi;
  #else
  const real prefact = 4.0*boil::pi;
  #endif
#endif

  const real volf_in  = u.bnd_flow( BndType::inlet() )
                      + u.bnd_flow( BndType::insert() )
                      + v_phase_change;

  /*---------------------------------------------+
  |  special outlet condition for bubble growth  |
  +---------------------------------------------*/
#ifdef TWODXY
  real deltz=u.zc(Comp::k(),ek(Comp::k()))-u.zc(Comp::k(),sk(Comp::k()));
  //boil::oout<<"deltz= "<<deltz<<"\n";
#elif defined TWODXZ
  real delty=u.yc(Comp::j(),ej(Comp::j()))-u.yc(Comp::j(),sj(Comp::j()));
#endif

  for(int iloop=1; iloop<=2; iloop++){
    real sumflx,ratio1;
    if(iloop==1){
      sumflx=0.0;
    } else {
      boil::cart.sum_real(&sumflx);
      //boil::oout<<"scale_out "<<volf_in<<" "<<sumflx<<"\n";
      if(sumflx==0.0){
        ratio1=1.0;
      } else {
        ratio1=volf_in/sumflx;
      }
    }
    for_m(m){
      for( int b=0; b<u.bc(m).count(); b++ ) {
      if( u.bc(m).type_decomp(b) ) continue;
        if( u.bc(m).type(b) == BndType::outlet() ) {
          Dir d = u.bc(m).direction(b);
          /* imin */
          if( m == Comp::u() && d == Dir::imin() ){
            for_vjk(u.bc(m).at(b),j,k){
              if( dom->ibody().off(si(m),j,k)) continue;
              real xx = u.xc(m,si(m)-1);
              real yy = u.yc(m,j);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(prefact*rr*deltz);
#elif defined TWODXZ
              real rr = sqrt(xx*xx + zz*zz);
              real uabs = volf_in/(prefact*rr*delty);
#elif defined AXISYM
              real rr = sqrt(xx*xx + zz*zz);           
              real uabs = volf_in/(prefact*rr*rr);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#endif
              if(iloop==1){
                sumflx += -dSx(m,si(m)-1,j,k) * uabs * xx / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][si(m)-1][j][k] = ratio1 * uabs * xx / rr;
                for(int ii=1; ii<=boil::BW; ii++)
                  u[m][si(m)-ii][j][k] = ratio1 * uabs * xx / rr;
              }
            }
          }
          /* imax */
          if( m == Comp::u() && d == Dir::imax() ){
            for_vjk(u.bc(m).at(b),j,k){
              if( dom->ibody().off(ei(m),j,k)) continue;
              real xx = u.xc(m,ei(m)+1);
              real yy = u.yc(m,j);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(prefact*rr*deltz);
#elif defined TWODXZ
              real rr = sqrt(xx*xx + zz*zz);
              real uabs = volf_in/(prefact*rr*delty);
#elif defined AXISYM
              real rr = sqrt(xx*xx + zz*zz);           
              real uabs = volf_in/(prefact*rr*rr);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#endif
              if(iloop==1){
                sumflx += dSx(m,ei(m)+1,j,k) * uabs * xx / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][ei(m)+1][j][k] = ratio1 * uabs * xx / rr;
                for(int ii=1; ii<=boil::BW; ii++)
                  u[m][ei(m)+ii][j][k] = ratio1 * uabs * xx / rr;
              }
            }
          }
          /* jmin */
          if( m == Comp::v() && d == Dir::jmin() ){
            for_vik(u.bc(m).at(b),i,k){
              if( dom->ibody().off(i,sj(m),k)) continue;
              real xx = u.xc(m,i);
              real yy = u.yc(m,sj(m)-1);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(prefact*rr*deltz);
              u[m][i][sj(m)-1][k] = uabs * yy / rr;
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
              u[m][i][sj(m)-1][k] = uabs * yy / rr;
#endif
#if (defined AXISYM || defined TWODXZ)
              boil::aout<<"momentum::scale_outlet_vel: not allowed!"
                        <<boil::endl;
              exit(0);
#endif
              if(iloop==1){
                sumflx += -dSy(m,i,sj(m)-1,k) * uabs * yy / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][i][sj(m)-1][k] = ratio1 * uabs * yy / rr;
                for(int jj=1; jj<=boil::BW; jj++)
                  u[m][i][sj(m)-jj][k] = ratio1 * uabs * yy / rr;
              }
            }
          }
          /* jmax */
          if( m == Comp::v() && d == Dir::jmax() ){
            for_vik(u.bc(m).at(b),i,k){
              if( dom->ibody().off(i,ej(m),k)) continue;
              real xx = u.xc(m,i);
              real yy = u.yc(m,ej(m)+1);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(prefact*rr*deltz);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#endif
#if (defined AXISYM || defined TWODXZ)
              boil::aout<<"momentum::scale_outlet_vel: not allowed!"
                        <<boil::endl;
              exit(0);
#endif
              if(iloop==1){
                sumflx += dSy(m,i,ej(m)+1,k) * uabs * yy / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][i][ej(m)+1][k] = ratio1 * uabs * yy / rr;
                for(int jj=1; jj<=boil::BW; jj++)
                  u[m][i][ej(m)+jj][k] = ratio1 * uabs * yy / rr;
              }
            }
          }
          /* kmin */
          if( m == Comp::w() && d == Dir::kmin() ){
            for_vij(u.bc(m).at(b),i,j){
              if( dom->ibody().off(i,j,sk(m))) continue;
              real xx = u.xc(m,i);
              real yy = u.yc(m,j);
              real zz = u.zc(m,sk(m)-1);
#ifdef TWODXY
              real uabs,rr;
              boil::aout<<"momentum::scale_outlet_vel: not allowed!"
                        <<boil::endl;
              exit(0);
#elif defined TWODXZ
              real rr = sqrt(xx*xx + zz*zz);
              real uabs = volf_in/(prefact*rr*delty);
#elif defined AXISYM
              real rr = sqrt(xx*xx + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#endif
              if(iloop==1){
                sumflx += -dSz(m,i,j,sk(m)-1) * uabs * zz / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][i][j][sk(m)-1] = ratio1 * uabs * zz / rr;
                for(int kk=1; kk<=boil::BW; kk++)
                  u[m][i][j][sk(m)-kk] = ratio1 * uabs * zz / rr;
              }
            }
          }
          /* kmax */
          if( m == Comp::w() && d == Dir::kmax() ){
            for_vij(u.bc(m).at(b),i,j){
              if( dom->ibody().off(i,j,ek(m))) continue;
              real xx = u.xc(m,i);
              real yy = u.yc(m,j);
              real zz = u.zc(m,ek(m)+1);
#ifdef TWODXY
              real uabs,rr;
              boil::aout<<"momentum::scale_outlet_vel: not allowed!"
                        <<boil::endl;
              exit(0);
#elif defined TWODXZ
              real rr = sqrt(xx*xx + zz*zz);
              real uabs = volf_in/(prefact*rr*delty);
#elif defined AXISYM
              real rr = sqrt(xx*xx + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(prefact*rr*rr);
#endif
              if(iloop==1){
                sumflx += dSz(m,i,j,ek(m)+1) * uabs * zz / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][i][j][ek(m)+1] = ratio1 * uabs * zz / rr;
                for(int kk=1; kk<=boil::BW; kk++)
                  u[m][i][j][ek(m)+kk] = ratio1 * uabs * zz / rr;
              }
            }
          }
        } /* bc == outlet */
      } /* bcs */
    } /* components */
  } /* loops */

  return;
}
