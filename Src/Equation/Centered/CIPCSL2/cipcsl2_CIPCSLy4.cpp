#include "cipcsl2.h"
#include <cmath>

/******************************************************************************/
void CIPCSL2::CIPCSLy4(const Vector & f, const Scalar & sy) {
/***************************************************************************//**
*  \brief Advance color-function in y-direction.      
*         j-face to cell.
*******************************************************************************/

  /* Define advection velocity: j-face */
  Comp m=Comp::j();
  for(int i=u->si(m)-1; i<=u->ei(m)+1; i++) {
  for(int j=u->sj(m)-1; j<=u->ej(m)+1; j++) {
  for(int k=u->sk(m)-1; k<=u->ek(m)+1; k++) {
    vel[i][j][k]=(*u)[m][i][j][k];
  }}}

  /* Reset delrho */
  for_aijk(i,j,k)
    delrho[i][j][k]=0.0;

  /* CIPCSL 1D */
  const real dt= time->dt();
  real y, den, d, dy, jsgn;
  real outd,outf;
  int jup,ist,ied,jst,jed,kst,ked;

  ist=si(); ied=ei(); jst=sj()+1; jed=ej(); kst=sk(); ked=ek();
  if( f.bc(m).type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0)                      jst=sj();
  if( f.bc(m).type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) != dom->dim(Comp::j())-1) jed=ej()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    y = -vel[i][j][k]*dt;
    jsgn = -copysign(1.0, y);
    jup = j - (int) jsgn;
    den = sy[i][j-1][k];
    dy = clr.dyc(j-1);
#ifdef DEBUG
    if(dy==0){
      std::cout<<"dy=0\n";
      exit(0);
    }
#endif
    if(y>0.0){
      den = sy[i][j][k];
      dy = clr.dyc(j);
    }
    d=-jsgn*dy ;
    funclxyz(outd,outf,f[m][i][j][k],f[m][i][jup][k],d,den,y,jsgn);
    fn[i][j][k]=outf;
    delrho[i][j][k]=outd;
  }}}

  /* b.c. for delrho */
  int jof,jof2;
  for( int b=0; b<clr.bc().count(); b++ ) {
    if ( clr.bc().type_decomp(b) ) continue;
    if ( clr.bc().direction(b)==Dir::jmin()
      || clr.bc().direction(b)==Dir::jmax() ){

      /* inlet or dirichlet */
      if ( clr.bc().type(b)==BndType::inlet()
        || clr.bc().type(b)==BndType::dirichlet()
        || clr.bc().type(b)==BndType::insert()){
        if (clr.bc().direction(b)==Dir::jmin()){
          jof=1; jof2=1;
        }
        if (clr.bc().direction(b)==Dir::jmax()){
          jof=0; jof2=0;
        }
        for_vijk( clr.bc().at(b), i,j,k ) {
          real dx=clr.dxc(i);
          real dz=clr.dzc(k);
          delrho[i][j+jof2][k]=phi[i][j][k]*vel[i][j+jof][k]*dt*dx*dz;
        }
      }

      /* outlet */
      if ( clr.bc().type(b)==BndType::outlet() ){
        if (clr.bc().direction(b)==Dir::jmin()){
          jof=1; jof2=1;
        }
        if (clr.bc().direction(b)==Dir::jmax()){
          jof=-1; jof2=0;
        }
        for_vijk( clr.bc().at(b), i,j,k ) {
          real dy=clr.dyc(j+jof);
          delrho[i][j+jof2][k]=clr[i][j+jof][k]/dy*vel[i][j+jof2][k]*dt;
          sum_outlet += -real(jof)*delrho[i][j+jof2][k];
          sum_outletm+= -real(jof)*(dV(i,j+jof,k)-clr[i][j+jof][k])
		         /dy*vel[i][j+jof2][k]*dt;
        }
      }
    }
  }

  /* update f */
  ist=si(); ied=ei(); jst=sj()+1; jed=ej(); kst=sk(); ked=ek();
  if( f.bc(m).type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0)                      jst=sj();
  if( f.bc(m).type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) != dom->dim(Comp::j())-1) jed=ej()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    real dy2=clr.dyc(j)+clr.dyc(j-1);
    f[m][i][j][k]=fn[i][j][k] 
              -fn[i][j][k]*(vel[i][j+1][k]-vel[i][j-1][k])*dt/dy2;
  }}}

  /* update sy */
  ist=si(); ied=ei(); jst=sj(); jed=ej(); kst=sk(); ked=ek();
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    sy[i][j][k]=sy[i][j][k]+(delrho[i][j][k]-delrho[i][j+1][k]);
  }}}
  return;
}
