#include "momentum.h"
//#define TWODXY

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
                u[m][si(m)-1][j][k] = ratio1 * uabs * xx / rr;
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
                u[m][ei(m)+1][j][k] = ratio1 * uabs * xx / rr;
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
                u[m][i][sj(m)-1][k] = ratio1 * uabs * yy / rr;
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
                u[m][i][ej(m)+1][k] = ratio1 * uabs * yy / rr;
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
                u[m][i][j][sk(m)-1] = ratio1 * uabs * zz / rr;
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
                u[m][i][j][ek(m)+1] = ratio1 * uabs * zz / rr;
              }
#endif
            }
         }
        }
      }
    }
  }
#endif

#if 0
  /*------------------------+
  |  nothing at the outlet  |
  +------------------------*/
  if( ratio > 2.0 || ratio <0.5 ) {
    /* bad, code duplication, this has been coppied from insert_bc */
    for_m(m){
    for( int b=0; b<u.bc(m).count(); b++ ) {

      if( u.bc(m).type_decomp(b) ) continue;

      if( u.bc(m).type(b) == BndType::outlet() ) {

        Dir d = u.bc(m).direction(b);

        if( m == Comp::u() && d == Dir::imin() ){
          //std::cout<<"dir= "<<d<<"m= "<<m<<"\n";
          for_vjk(u.bc(m).at(b),j,k) u[m][si(m)-1][j][k] = -volf_in/asum;
        }
        if( m == Comp::u() && d == Dir::imax() ){
          //std::cout<<"dir= "<<d<<"m= "<<m<<"\n";
          for_vjk(u.bc(m).at(b),j,k) u[m][ei(m)+1][j][k] = volf_in/asum;
        }

        if( m == Comp::v() && d == Dir::jmin() ){
          //std::cout<<"dir= "<<d<<"m= "<<m<<"\n";
          for_vik(u.bc(m).at(b),i,k) u[m][i][sj(m)-1][k] = -volf_in/asum;
        }
        if( m == Comp::v() && d == Dir::jmax() ){
          //std::cout<<"dir= "<<d<<"m= "<<m<<"\n";
          for_vik(u.bc(m).at(b),i,k) u[m][i][ej(m)+1][k] = volf_in/asum;
        }

        if( m == Comp::w() && d == Dir::kmin() ){
          //std::cout<<"dir= "<<d<<"m= "<<m<<"\n";
          for_vij(u.bc(m).at(b),i,j) u[m][i][j][sk(m)-1] = -volf_in/asum;
        }
        if( m == Comp::w() && d == Dir::kmax() ){
          //std::cout<<"dir= "<<d<<"m= "<<m<<"\n";
          for_vij(u.bc(m).at(b),i,j) u[m][i][j][ek(m)+1] = volf_in/asum;
        }
      }
    } /* b */
    }
  }
  /*--------------------------+
  |  something at the outlet  |
  +--------------------------*/
  else {
    for_m(m) {              
      /* again code duplication, this has been coppied from insert_bc */
      for( int b=0; b<u.bc(m).count(); b++ ) {

        if( u.bc(m).type_decomp(b) ) continue;

        if( u.bc(m).type(b) == BndType::outlet() ) {

          Dir d = u.bc(m).direction(b);

          if( d == Dir::imin() ) 
            for_vjk(u.bc(m).at(b),j,k) u[m][si(m)-1][j][k] *= ratio;    
          if( d == Dir::imax() ) 
            for_vjk(u.bc(m).at(b),j,k) u[m][ei(m)+1][j][k] *= ratio;    

          if( d == Dir::jmin() ) 
            for_vik(u.bc(m).at(b),i,k) u[m][i][sj(m)-1][k] *= ratio;    
          if( d == Dir::jmax() ) 
            for_vik(u.bc(m).at(b),i,k) u[m][i][ej(m)+1][k] *= ratio;    

          if( d == Dir::kmin() ) 
            for_vij(u.bc(m).at(b),i,j) u[m][i][j][sk(m)-1] *= ratio;    
            //for_vij(u.bc(m).at(b),i,j) u[m][i][j][sk(m)+1] *= ratio;    
          if( d == Dir::kmax() ) 
            for_vij(u.bc(m).at(b),i,j) u[m][i][j][ek(m)+1] *= ratio;    
            //for_vij(u.bc(m).at(b),i,j) u[m][i][j][ek(m)-1] *= ratio;    

        }
      } /* b */
    } /* m */
  }
#endif
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_scale_out.cpp,v 1.1 2015/02/19 09:41:42 sato Exp $'/
+-----------------------------------------------------------------------------*/
