#include "cipcsl2.h"
#include <iomanip>
#define WALL
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::smear(Scalar & sca) {
/***************************************************************************//**
* \brief  Sharpen color function by modified Olsson's method.
*            input & output: sca
*            temporary: stmp
*******************************************************************************/
  real dsca_max=0.0;  // l_infinity for iterative cal.
  const real epss = dxmin*dxmin;
  const real dtau=0.125;

  /* flag of boundary conditions for stmp */
  bool imin, imax, jmin, jmax, kmin, kmax;
  imin = imax = jmin = jmax = kmin = kmax = false;
  for( int b=0; b<phi.bc().count(); b++ ) {
    if( phi.bc().type_decomp(b) ) continue;
    if( phi.bc().type(b)==BndType::dirichlet() ||
        phi.bc().type(b)==BndType::inlet()     ||
        phi.bc().type(b)==BndType::neumann()   ||
        phi.bc().type(b)==BndType::wall()      ||
        phi.bc().type(b)==BndType::outlet()    ||
        phi.bc().type(b)==BndType::symmetry()  ||
        phi.bc().type(b)==BndType::insert()    ) {
      if(phi.bc().direction(b)==Dir::imin()) imin=true;
      if(phi.bc().direction(b)==Dir::imax()) imax=true;
      if(phi.bc().direction(b)==Dir::jmin()) jmin=true;
      if(phi.bc().direction(b)==Dir::jmax()) jmax=true;
      if(phi.bc().direction(b)==Dir::kmin()) kmin=true;
      if(phi.bc().direction(b)==Dir::kmax()) kmax=true;
    }
  }

#ifdef WALL
  /* flag of boundary conditions for stmp & wall */
  bool iminw, imaxw, jminw, jmaxw, kminw, kmaxw;
  iminw = imaxw = jminw = jmaxw = kminw = kmaxw = false;
  for( int b=0; b<phi.bc().count(); b++ ) {
    if( phi.bc().type_decomp(b) ) continue;
    if( phi.bc().type(b)==BndType::wall() ){
      if(phi.bc().direction(b)==Dir::imin()) iminw=true;
      if(phi.bc().direction(b)==Dir::imax()) imaxw=true;
      if(phi.bc().direction(b)==Dir::jmin()) jminw=true;
      if(phi.bc().direction(b)==Dir::jmax()) jmaxw=true;
      if(phi.bc().direction(b)==Dir::kmin()) kminw=true;
      if(phi.bc().direction(b)==Dir::kmax()) kmaxw=true;
    }
  }
#endif

  /* iterate */
  for(int it=0; it<itsmear; it++) {

#ifdef DEBUG
  boil::oout<<"it="<<it<<" "<<itsmear<<"\n";
#endif

    stmp = 0.0;

    // x-direction
    //for_vmijk((*u),m,i,j,k){  //don't use vmijk. wall will be skipped!
    for(int i=si(); i<=ei()+1; i++)
      for(int j=sj(); j<=ej(); j++)
        for(int k=sk(); k<=ek(); k++){
          if (dom->ibody().off(Comp::u(),i,j,k)) continue;
          real fluxd, ldt;
          // diffusion
          fluxd = epss * (sca[i][j][k]-sca[i-1][j][k])*dSx(i,j,k)/sca.dxw(i);
          ldt  = dtau;

          if(i==si()   && imin) ldt=0.0;
          if(i==ei()+1 && imax) ldt=0.0;
          if(i==si()+1 && iminw) ldt=0.0;
          if(i==ei()   && imaxw) ldt=0.0;
          if(dom->ibody().off(Comp::u(),i-1,j,k)) ldt=0.0; // adjascent wall
          if(dom->ibody().off(Comp::u(),i+1,j,k)) ldt=0.0; // adjascent wall

          // add flux
          stmp[i-1][j][k] -= ( - fluxd)*ldt;
          stmp[i  ][j][k] += ( - fluxd)*ldt;
        }

    // y-direction
    //for_vmijk((*u),m,i,j,k){  //don't use vmijk. wall will be skipped!
    for(int i=si(); i<=ei(); i++)
      for(int j=sj(); j<=ej()+1; j++)
        for(int k=sk(); k<=ek(); k++){
          if (dom->ibody().off(Comp::v(),i,j,k)) continue;
          real fluxd, ldt;
          // diffusion
          fluxd = epss * (sca[i][j][k]-sca[i][j-1][k])*dSy(i,j,k)/sca.dys(j);
          ldt  = dtau;

          if(j==sj()   && jmin) ldt=0.0;
          if(j==ej()+1 && jmax) ldt=0.0;
          if(j==sj()+1 && jminw) ldt=0.0;
          if(j==ej()   && jmaxw) ldt=0.0;
          if(dom->ibody().off(Comp::v(),i,j-1,k)) ldt=0.0; // adjascent wall
          if(dom->ibody().off(Comp::v(),i,j+1,k)) ldt=0.0; // adjascent wall

          // add flux
          stmp[i][j-1][k] -= ( - fluxd)*ldt;
          stmp[i][j  ][k] += ( - fluxd)*ldt;
        }

    // z-direction
    //for_vmijk((*u),m,i,j,k){  //don't use vmijk. wall will be skipped!
    for(int i=si(); i<=ei(); i++)
      for(int j=sj(); j<=ej(); j++)
        for(int k=sk(); k<=ek()+1; k++){
          if (dom->ibody().off(Comp::w(),i,j,k)) continue;
          real fluxd, ldt;
          // diffusion
          fluxd = epss * (sca[i][j][k]-sca[i][j][k-1])*dSz(i,j,k)/sca.dzb(k);
          ldt  = dtau;

          if(k==sk()   && kmin) ldt=0.0;
          if(k==ek()+1 && kmax) ldt=0.0;
          if(k==sk()+1 && kminw) ldt=0.0;
          if(k==ek()   && kmaxw) ldt=0.0;
          if(dom->ibody().off(Comp::w(),i,j,k-1)) ldt=0.0; // adjascent wall
          if(dom->ibody().off(Comp::w(),i,j,k+1)) ldt=0.0; // adjascent wall

          // add flux
          stmp[i][j][k-1] -= ( - fluxd)*ldt;
          stmp[i][j][k  ] += ( - fluxd)*ldt;
        }

    /* advance */
    for_ijk(i,j,k){
      real dsca = 1.0/dV(i,j,k) * stmp[i][j][k];
      sca[i][j][k] += dsca;
    }

    /* boundary condition */
    //sca.bnd_update();
    bdcond(sca);
    sca.exchange_all();

    /*--------------------+
    |  immersed boundary  |
    +--------------------*/
    ib_ext_scalar(sca);
    ib_bdcond(sca); 

  }

  return;
}
