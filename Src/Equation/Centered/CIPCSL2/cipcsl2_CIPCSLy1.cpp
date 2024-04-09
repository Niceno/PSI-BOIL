#include "cipcsl2.h"
#include <cmath>

/******************************************************************************/
void CIPCSL2::CIPCSLy1(const Scalar & f, const Scalar & sy) {
/***************************************************************************//**
*  \brief Advance color-function in y-direction.      
*         node to j-edge.
*******************************************************************************/

   /* Define advection velocity: node */
   Comp m=Comp::v();
   for(int i=u->si(m);   i<=u->ei(m)+1; i++) {
   for(int j=u->sj(m)-1; j<=u->ej(m)+1; j++) {
   for(int k=u->sk(m);   k<=u->ek(m)+1; k++) {
     vel[i][j][k]=((*u)[m][i  ][j][k]+(*u)[m][i][j][k-1]+
                   (*u)[m][i-1][j][k]+(*u)[m][i][j][k-1])/4.0;
   }}}

  /* reset delrho */
  for_aijk(i,j,k)
    delrho[i][j][k]=0.0;

  /* CIPCSL 1D */
  const real dt= time->dt();
  real y, den, d, dy, jsgn;
  real outd,outf;
  int jup,ist,ied,jst,jed,kst,ked;

  ist=si()+1; ied=ei(); jst=sj()+1; jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0)                      ist=si();
  if( f.bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) != dom->dim(Comp::i())-1) ied=ei()+1;
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
    y=-vel[i][j][k]*dt;
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
    funclxyz(outd,outf,f[i][j][k],f[i][jup][k],d,den,y,jsgn);
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
          jof=1;
          jof2=1;
        }
        if (clr.bc().direction(b)==Dir::jmax()){
          jof=0;
          jof2=0;
        }
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()+1; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()  ; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()+1; k++){
          real vbnd=0.25*(phi[i][j][k  ]+phi[i-1][j][k]
                         +phi[i][j][k-1]+phi[i-1][j][k-1]);
          delrho[i][j+jof2][k]=vbnd*vel[i][j+jof][k]*dt;
        }}}
      }

      /* outlet */
      if ( clr.bc().type(b)==BndType::outlet() ){
        if (clr.bc().direction(b)==Dir::jmin()){
          jof=1;
          jof2=1;
        }
        if (clr.bc().direction(b)==Dir::jmax()){
          jof=-1;
          jof2=0;
        }
        for(int i=clr.bc().at(b).si(); i<=clr.bc().at(b).ei()+1; i++){
        for(int j=clr.bc().at(b).sj(); j<=clr.bc().at(b).ej()  ; j++){
        for(int k=clr.bc().at(b).sk(); k<=clr.bc().at(b).ek()+1; k++){
          real dy=clr.dyc(j+jof);
          delrho[i][j+jof2][k]=sy[i][j+jof][k]/dy*vel[i][j+jof2][k]*dt;
        }}}
      }
    }
  }

  /* update f */
  ist=si()+1; ied=ei(); jst=sj()+1; jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0)                      ist=si();
  if( f.bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) != dom->dim(Comp::i())-1) ied=ei()+1;
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
    real dy2=clr.dyc(j)+clr.dyc(j-1);
    f[i][j][k]=fn[i][j][k] 
              -fn[i][j][k]*(vel[i][j+1][k]-vel[i][j-1][k])*dt/dy2;
  }}}

  /* update sy */
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
    sy[i][j][k]=sy[i][j][k]+(delrho[i][j][k]-delrho[i][j+1][k]) ;
  }}}
  return;
}
