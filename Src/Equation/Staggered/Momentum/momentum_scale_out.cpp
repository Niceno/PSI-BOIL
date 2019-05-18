#include "momentum.h"
#define TWODXY

/******************************************************************************/
void Momentum::scale_out() {
/*----------------------------------+ 
|  set bulk velocity at the outlet  |
+----------------------------------*/

  Comp m;

  /*---------------------------------+
  |  if outlet does not exist, exit  |
  +---------------------------------*/
  if( 
      !u.bc( Comp::u() ).exists( BndType::outlet() ) &&
      !u.bc( Comp::v() ).exists( BndType::outlet() ) && 
      !u.bc( Comp::w() ).exists( BndType::outlet() ) 
    ) return;

  /*----------------------------+
  |  get volume flux the inlet  |
  +----------------------------*/
  // v_phase_change() should be called from main.cpp
  const real volf_in  = volf_bct( BndType::inlet() ) + v_phase_change;
  OPR(volf_in);
  //assert(volf_in != 0.0);

  /*----------------------------------------+
  |  get current volume flux at the outlet  |
  +----------------------------------------*/
  real ax=0.0, ay=0.0, az=0.0;
  const real volf_out = volf_bct( BndType::outlet(), &ax, &ay, &az );
  //real am = boil::maxr(ax, ay, az);
  //std::cout<<"ax,ay,az= "<<ax<<" "<<ay<<" "<<" "<<az<<"\n";
  real asum = ax + ay + az;

  real ratio;
  if(volf_out==0.0 && volf_in==0.0){
    ratio=0.0;
  } else if(volf_out==0.0 && volf_in!=0.0){
    ratio=1.0e+12;
  } else {
    ratio= -volf_in / volf_out;
  }

#if 1
  /*---------------------------------------------+
  |  special outlet condition for bubble growth  |
  +---------------------------------------------*/
#ifdef TWODXY
  real deltz=u.zc(Comp::k(),ek(Comp::k()))-u.zc(Comp::k(),sk(Comp::k()));
  //boil::oout<<"deltz= "<<deltz<<"\n";
#endif
  real pi=acos(-1.0);
  //boil::oout<<"pi= "<<pi<<"\n";

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
              real xx = u.xc(m,si(m)-1);
              real yy = u.yc(m,j);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(2.0*pi*rr*deltz);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(4.0*pi*rr*rr);
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
              real xx = u.xc(m,ei(m)+1);
              real yy = u.yc(m,j);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(2.0*pi*rr*deltz);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(4.0*pi*rr*rr);
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
              real xx = u.xc(m,i);
              real yy = u.yc(m,sj(m)-1);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(2.0*pi*rr*deltz);
              u[m][i][sj(m)-1][k] = uabs * yy / rr;
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(4.0*pi*rr*rr);
              u[m][i][sj(m)-1][k] = uabs * yy / rr;
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
              real xx = u.xc(m,i);
              real yy = u.yc(m,ej(m)+1);
              real zz = u.zc(m,k);
#ifdef TWODXY
              real rr = sqrt(xx*xx + yy*yy);
              real uabs = volf_in/(2.0*pi*rr*deltz);
#else
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(4.0*pi*rr*rr);
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
              real xx = u.xc(m,i);
              real yy = u.yc(m,j);
              real zz = u.zc(m,sk(m)-1);
#ifndef TWODXY
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(4.0*pi*rr*rr);
              if(iloop==1){
                sumflx += -dSz(m,i,j,sk(m)-1) * uabs * zz / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][i][j][sk(m)-1] = ratio1 * uabs * zz / rr;
                for(int kk=1; kk<=boil::BW; kk++)
                  u[m][i][j][sk(m)-kk] = ratio1 * uabs * zz / rr;
              }
#endif
            }
          }
          /* kmax */
          if( m == Comp::w() && d == Dir::kmax() ){
            for_vij(u.bc(m).at(b),i,j){
              real xx = u.xc(m,i);
              real yy = u.yc(m,j);
              real zz = u.zc(m,ek(m)+1);
#ifndef TWODXY
              real rr = sqrt(xx*xx + yy*yy + zz*zz);
              real uabs = volf_in/(4.0*pi*rr*rr);
              if(iloop==1){
                sumflx += dSz(m,i,j,ek(m)+1) * uabs * zz / rr;
              } else {
                /* value is only copied. The solid angle expansion
                   in the buffers should be also considered...*/
                //u[m][i][j][ek(m)+1] = ratio1 * uabs * zz / rr;
                for(int jj=1; jj<=boil::BW; jj++)
                  u[m][i][j][ek(m)+kk] = ratio1 * uabs * zz / rr;
              }
#endif
            }
         }
        }
      }
    }
  }
#endif
}
