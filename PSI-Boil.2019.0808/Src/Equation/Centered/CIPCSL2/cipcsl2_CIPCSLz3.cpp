#include "cipcsl2.h"
#include <cmath>

/******************************************************************************/
void CIPCSL2::CIPCSLz3(const Scalar & f, const Vector & sz) {
/***************************************************************************//**
*  \brief Advance color-function in z-direction.      
*         j-edge to i-face.
*******************************************************************************/

  /* Define advection velocity: j-line */
  Comp m=Comp::w();
  Comp mc=Comp::u();
  for(int i=u->si(m);   i<=u->ei(m)+1; i++) {
  for(int j=u->sj(m)-1; j<=u->ej(m)+1; j++) {
  for(int k=u->sk(m)-1; k<=u->ek(m)+1; k++) {
    vel[i][j][k]=((*u)[m][i][j][k]+(*u)[m][i-1][j][k])/2.0;
  }}}
  /* EXTENDED BUFFERS HINT:
     Lines 14-16 could probably be replaced with this:
  for(int i=u->si(m);   i<=u->ei(m)+1; i++) {
  for(int j=u->sj(m)-1; j<=u->ej(m)+1; j++) {
  for(int k=u->sk(m)-1; k<=u->ek(m)+1; k++) {
  */

  /* Reset delrho */
  for_aijk(i,j,k)
    delrho[i][j][k]=0.0;
  /* EXTENDED BUFFERS HINT:
     Lines 27-29 could probably be replaced with this:
  for_aijk(i,j,k)
    delrho[i][j][k]=0.0;
  */

  /* CIPCSL 1D */
  const real dt= time->dt();
  real z, den, d, dz, ksgn;
  real outd,outf;
  int kup,ist,ied,jst,jed,kst,ked;

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
    z = -vel[i][j][k]*dt;
    ksgn = -copysign(1.0, z);
    kup = k - (int) ksgn;
    den = sz[mc][i][j][k-1];
    dz = clr.dzc(k-1);
#ifdef DEBUG
    if(dz==0){
      std::cout<<"CIPCSLz3:dz=0\n";
      exit(0);
    }
#endif
    if(z>0.0){
      den = sz[mc][i][j][k];
      dz = clr.dzc(k);
    }
    d=-ksgn*dz;
    funclxyz(outd,outf,f[i][j][k],f[i][j][kup],d,den,z,ksgn);
    fn[i][j][k]=outf;
    delrho[i][j][k]=outd;
  }}}

  /* b.c. for delrho */
  int kof,kof2;
  for( int b=0; b<clr.bc().count(); b++ ) {
    if ( clr.bc().type_decomp(b) ) continue;
    if ( clr.bc().direction(b)==Dir::kmin()
      || clr.bc().direction(b)==Dir::kmax() ){
      /* inlet or dirichlet */
      if ( clr.bc().type(b)==BndType::inlet()
        || clr.bc().type(b)==BndType::dirichlet()
        || clr.bc().type(b)==BndType::insert()){
        if (clr.bc().direction(b)==Dir::kmin()){
          kof=1; kof2=1;
        }
        if (clr.bc().direction(b)==Dir::kmax()){
          kof=0; kof2=0;
        }
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()+1; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()  ; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()  ; k++){
          real dy=clr.dyc(j);
          real vbnd=0.5*(phi[i][j][k]+phi[i-1][j][k]);
          delrho[i][j][k+kof2]=vbnd*dy*vel[i][j][k+kof]*dt;
        }}}
      }
      /* outlet */
      if ( clr.bc().type(b)==BndType::outlet() ){
        if (clr.bc().direction(b)==Dir::kmin()){
          kof=1; kof2=1;
        }
        if (clr.bc().direction(b)==Dir::kmax()){
          kof=-1; kof2=0;
        }
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()+1; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()  ; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()  ; k++){
          real dz=clr.dzc(k+kof);
          delrho[i][j][k+kof2]=sz[mc][i][j][k+kof]/dz*vel[i][j][k+kof2]*dt;
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
    real dz2=clr.dzc(k)+clr.dzc(k-1);
    f[i][j][k]=fn[i][j][k]
              -fn[i][j][k]*(vel[i][j][k+1]-vel[i][j][k-1])*dt/dz2;
  }}}

  /* update sz */
  ist=si()+1; ied=ei(); jst=sj(); jed=ej(); kst=sk(); ked=ek();
  if( f.bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0)                      ist=si();
  if( f.bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) != dom->dim(Comp::i())-1) ied=ei()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){    
    sz[mc][i][j][k]=sz[mc][i][j][k]+(delrho[i][j][k]-delrho[i][j][k+1]) ;
  }}}
  return;
}
