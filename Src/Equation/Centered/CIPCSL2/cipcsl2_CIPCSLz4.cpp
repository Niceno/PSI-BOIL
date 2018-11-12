#include "cipcsl2.h"
#include <cmath>
using namespace std;

/**************** Advection in Z-direction ************************************/
void CIPCSL2::CIPCSLz4(const Vector & f, const Scalar & sz) {
/***************************************************************************//**
*  \brief Advance color-function in z-direction.
*         k-face to cell.
*******************************************************************************/

  /* Define advection velocity: k-face */
  Comp m=Comp::w();
  for(int i=0; i<=u->ni(m)-1; i++) {
  for(int j=0; j<=u->nj(m)-1; j++) {
  for(int k=0; k<=u->nk(m)-1; k++) {
    vel[i][j][k]=(*u)[m][i][j][k];
  }}}

  /* Reset delrho */
  for(int i=0; i<=ni(); i++)
  for(int j=0; j<=nj(); j++)
  for(int k=0; k<=nk(); k++)
    delrho[i][j][k]=0.0;

  /* CIPCSL 1D */
  const real dt= time->dt();
  real z, den, d, dz, ksgn;
  real outd,outf;
  int kup,ist,ied,jst,jed,kst,ked;

  ist=si(); ied=ei(); jst=sj(); jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc(m).type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk();
  if( f.bc(m).type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    z = -vel[i][j][k]*dt;
    ksgn = -copysign(1.0, z);
    kup = k - (int) ksgn;
    den = sz[i][j][k-1];
    dz = clr.dzc(k-1);
#ifdef DEBUG
    if(dz==0){
      std::cout<<"dz=0\n";
      exit(0);
    }
#endif
    if(z>0.0){
      den = sz[i][j][k];
      dz = clr.dzc(k);
    }
    d=-ksgn*dz;
    funclxyz(outd,outf,f[m][i][j][k],f[m][i][j][kup],d,den,z,ksgn);
    fn[i][j][k]=outf;
    delrho[i][j][k]=outd;
#if 0
    if(i==36&&j==29&&(k==2||k==10)) {
      cout<<"loop: "<<k<<" "<<outd<<" "<<outf<<" "<<f[m][i][j][k]<<" "<<f[m][i][j][kup]<<" "<<d<<" "<<den<<" "<<z<<" "<<ksgn<<"\n";
    }
#endif
  }}}
#if 0
  cout<<"Lz4:delrho "<<delrho[36][29][2]<<" "<<delrho[36][29][10]<<"\n";
#endif

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
        for_vijk( clr.bc().at(b), i,j,k ) {
          real dx=clr.dxc(i);
          real dy=clr.dyc(j);
          delrho[i][j][k+kof2]=phi[i][j][k]*vel[i][j][k+kof]*dt*dx*dy;
        }
      }

      /* outlet */
      if ( clr.bc().type(b)==BndType::outlet() ){
        if (clr.bc().direction(b)==Dir::kmin()){
          kof=1; kof2=1;
        }
        if (clr.bc().direction(b)==Dir::kmax()){
          kof=-1; kof2=0;
        }
        for_vijk( clr.bc().at(b), i,j,k ) {
          real dz=clr.dzc(k+kof);
          delrho[i][j][k+kof2]=clr[i][j][k+kof]/dz*vel[i][j][k+kof2]*dt;
          sum_outlet += -real(kof)*delrho[i][j][k+kof2];
          sum_outletm+= -real(kof)*(dV(i,j,k+kof)-clr[i][j][k+kof])
                        /dz*vel[i][j][k+kof2]*dt;
        }
      }
    }
  }

  /* update f */
  ist=si(); ied=ei(); jst=sj(); jed=ej(); kst=sk()+1; ked=ek();
  if( f.bc(m).type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk();
  if( f.bc(m).type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    real dz2=clr.dzc(k)+clr.dzc(k-1);
    f[m][i][j][k]=fn[i][j][k]
              -fn[i][j][k]*(vel[i][j][k+1]-vel[i][j][k-1])*dt/dz2;
  }}}

#if 0
  cout<<"Lz4: "<<delrho[36][29][1]<<" "<<delrho[36][29][2]<<" "<<vel[36][29][2]<<"\n";
  cout<<"Lz4: "<<delrho[36][29][9]<<" "<<delrho[36][29][10]<<" "<<vel[36][29][10]<<"\n";
#endif

  /* update sz */
  ist=si(); ied=ei(); jst=sj(); jed=ej(); kst=sk(); ked=ek();
  for(int i=ist; i<=ied; i++){
  for(int j=jst; j<=jed; j++){
  for(int k=kst; k<=ked; k++){
    sz[i][j][k]=sz[i][j][k]+(delrho[i][j][k]-delrho[i][j][k+1]);
  }}}
  return;
}
