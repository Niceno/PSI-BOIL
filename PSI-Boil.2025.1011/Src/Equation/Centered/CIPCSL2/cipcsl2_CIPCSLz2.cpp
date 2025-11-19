#include "cipcsl2.h"
#include <cmath>

/******************************************************************************/
void CIPCSL2::CIPCSLz2(const Scalar & f, const Vector & sz) {
/***************************************************************************//**
*  \brief Advance color-function in z-direction.
*         i-edge to j-face.
*******************************************************************************/

  /* Define advection velocity: i-line */
  Comp m=Comp::w();
  Comp mc=Comp::v();
  for(int i=u->si(m)-1; i<=u->ei(m)+1; i++) {
  for(int j=u->sj(m);   j<=u->ej(m)+1; j++) {
  for(int k=u->sk(m)-1; k<=u->ek(m)+1; k++) {
    vel[i][j][k]=((*u)[m][i][j][k]+(*u)[m][i][j-1][k])/2.0;
  }}}

  /* Reset delrho */
  for_aijk(i,j,k)
    delrho[i][j][k]=0.0;

  /* CIPCSL 1D */
  const real dt= time->dt();
  real z, den, d, dz, ksgn;
  real outd,outf;
  int kup,ist,ied,jst,jed,kst,ked;

  ist=si(); ied=ei(); jst=sj()+1; jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0)                      jst=sj();
  if( f.bc().type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) != dom->dim(Comp::j())-1) jed=ej()+1;
  if( f.bc().type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk();
  if( f.bc().type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    dz=clr.dzc(k);
    z = -vel[i][j][k]*dt;
    ksgn = -copysign(1.0, z);
    kup = k - (int) ksgn;
    den = sz[mc][i][j][k-1];
    dz = clr.dzc(k-1);
#ifdef DEBUG
    if(dz==0){
      std::cout<<"CIPCSLLz2:dz=0\n";
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
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()  ; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()+1; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()  ; k++){
          real dx=clr.dxc(i);
          real vbnd=0.5*(phi[i][j][k]+phi[i][j-1][k]);
          delrho[i][j][k+kof2]=vbnd*dx*vel[i][j][k+kof]*dt;
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
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()  ; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()+1; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()  ; k++){
          real dz=clr.dzc(k+kof);
          delrho[i][j][k+kof2]=sz[mc][i][j][k+kof]/dz*vel[i][j][k+kof2]*dt;
        }}}
      }
    }
  }

  /* update f */
  ist=si(); ied=ei(); jst=sj()+1; jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0)                      jst=sj();
  if( f.bc().type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) != dom->dim(Comp::j())-1) jed=ej()+1;
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
  ist=si(); ied=ei(); jst=sj()+1; jed=ej(); kst=sk(); ked=ek();
  if( f.bc().type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0)                      jst=sj();
  if( f.bc().type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) != dom->dim(Comp::j())-1) jed=ej()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){    
    sz[mc][i][j][k]=sz[mc][i][j][k]+(delrho[i][j][k]-delrho[i][j][k+1]) ;
  }}}
  return;
}
