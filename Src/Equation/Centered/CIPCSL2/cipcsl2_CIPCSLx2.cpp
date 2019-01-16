#include "cipcsl2.h"
#include <cmath>

/******************************************************************************/
void CIPCSL2::CIPCSLx2(const Scalar & f, const Vector & sx) {
/***************************************************************************//**
*  \brief Advance color-function in x-direction.
*         j-edge to k-face.
*******************************************************************************/

  /* Define advection velocity: j-line */
  Comp m=Comp::u();
  Comp mc=Comp::w();
  for(int i=0; i<=u->ni(m)-1; i++) {
  for(int j=0; j<=u->nj(m)-1; j++) {
  for(int k=1; k<=u->nk(m)-1; k++) {
    vel[i][j][k]=((*u)[m][i][j][k]+(*u)[m][i][j][k-1])/2.0;
  }}}
  /* EXTENDED BUFFERS HINT:
     Lines 14-16 could probably be replaced with this:
  for(int i=u->si(m)-1; i<=u->ei(m)+1; i++) {
  for(int j=u->sj(m)-1; j<=u->ej(m)+1; j++) {
  for(int k=u->sk(m);   k<=u->ek(m)+1; k++) {
  */

  /* Reset delrho */
  for(int i=0; i<=ni(); i++)
  for(int j=0; j<=nj(); j++)
  for(int k=0; k<=nk(); k++)
    delrho[i][j][k]=0.0;
  /* EXTENDED BUFFERS HINT:
     Lines 27-29 could probably be replaced with this:
  for_aijk(i,j,k)
    delrho[i][j][k]=0.0;
  */

  /* CIPCSL 1D */
  const real dt= time->dt();
  real x, den, d, dx, isgn;
  real outd,outf;
  int iup,ist,ied,jst,jed,kst,ked;
  ist=si()+1; ied=ei(); jst=sj(); jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0)                      ist=si();
  if( f.bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) != dom->dim(Comp::i())-1) ied=ei()+1;
  if( f.bc().type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk();
  if( f.bc().type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    x=-vel[i][j][k]*dt;
    isgn = -copysign(1.0, x);
    iup = i - (int) isgn;
    den = sx[mc][i-1][j][k];
    dx = clr.dxc(i-1);
#ifdef DEBUG
    if(dx==0){
      std::cout<<"dx=0\n";
      exit(0);
    }
#endif
    if(x>0.0){
      den = sx[mc][i][j][k];
      dx = clr.dxc(i);
    }
    d=-isgn*dx;
    funclxyz(outd,outf,f[i][j][k],f[iup][j][k],d,den,x,isgn);
    fn[i][j][k]=outf;
    delrho[i][j][k]=outd;
  }}}

  /* b.c. for delrho */
  int iof,iof2;
  for( int b=0; b<clr.bc().count(); b++ ) {
    if ( clr.bc().type_decomp(b) ) continue;
    if ( clr.bc().direction(b)==Dir::imin()
      || clr.bc().direction(b)==Dir::imax() ){
      /* inlet or dirichlet */
      if ( clr.bc().type(b)==BndType::inlet()
        || clr.bc().type(b)==BndType::dirichlet()
        || clr.bc().type(b)==BndType::insert()){
        if (clr.bc().direction(b)==Dir::imin()){
          iof=1; iof2=1;
        }
        if (clr.bc().direction(b)==Dir::imax()){
          iof=0; iof2=0;
        }
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()  ; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()  ; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()+1; k++){
          real dy=clr.dyc(j);
          real vbnd=0.5*(phi[i][j][k  ]+phi[i][j][k-1]);
          delrho[i+iof2][j][k]=vbnd*dy*vel[i+iof][j][k]*dt;
        }}}
      }
      /* outlet */
      if ( clr.bc().type(b)==BndType::outlet() ){
        if (clr.bc().direction(b)==Dir::imin()){
          iof=1; iof2=1;
        }
        if (clr.bc().direction(b)==Dir::imax()){
          iof=-1; iof2=0;
        }
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()  ; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()  ; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()+1; k++){
          real dx=clr.dxc(i+iof);
          delrho[i+iof2][j][k]=sx[mc][i+iof][j][k]/dx*vel[i+iof2][j][k]*dt;
        }}}
      }
    }
  }

  /* update f */
  ist=si()+1; ied=ei(); jst=sj(); jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0)                      ist=si();
  if( f.bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) != dom->dim(Comp::i())-1) ied=ei()+1;
  if( f.bc().type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk();
  if( f.bc().type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    real dx2=clr.dxc(i)+clr.dxc(i-1);
    f[i][j][k]=fn[i][j][k] 
              -fn[i][j][k]*(vel[i+1][j][k]-vel[i-1][j][k])*dt/dx2;
  }}}

  /* update sx */
  ist=si(); ied=ei(); jst=sj(); jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk();
  if( f.bc().type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    dx=clr.dxc(i);
    sx[mc][i][j][k]=sx[mc][i][j][k]+(delrho[i][j][k]-delrho[i+1][j][k]);
  }}}
  return;
}
