#include "cipcsl2.h"
#include <iomanip>
#define WALL
//#define DEBUG
using namespace std;

//real beta(const real a1, const real a2, const bool local);

/******************************************************************************/
bool CIPCSL2::sharpen(Scalar & sca, const real epss, const int itmax, const bool local) {
/***************************************************************************//**
* \brief  Sharpen color function by modified Olsson's method.
*            input & output: sca
*            temporary: stmp
*******************************************************************************/
  real dsca_max=0.0;  // l_infinity for iterative cal.

  bool ibm=false;
  if(dom->ibody().nccells() > 0) {
    ibm=true;
  }

#ifdef DEBUG
  boil::oout<<"sharpen:ibm= "<<ibm<<"\n";
#endif

  //for_ijk(i,j,k){
  for_aijk(i,j,k){
    real dxtmp = max(sca.dxc(i),sca.dyc(j));
         dxtmp = max(sca.dzc(k),dxtmp);
    atmp[i][j][k] = dxtmp;
    if (dxtmp<=0.0) {
      fn[i][j][k] = 0.0;
    } else {
      /* diffusion number = eps * dt / (dx)^2 */
      fn[i][j][k] = min(0.5*time->dt(),0.125*dxtmp*dxtmp/(epss*dxtmp));
      //fn  [i][j][k] = 0.125*dxtmp*dxtmp/(epss*dxtmp);
    }
  }
  //insert_bc(fn);
  fn.exchange();
  //insert_bc(atmp);
  atmp.exchange();

#ifdef DEBUG
  boil::oout<<"sharpen:insert_bc \n";
#endif

#if 0
  if(time->current_step()==1) {
    boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
    boil::plot->plot(dist,nx,ny,nz, "dist-nx-ny-nz", time->current_step());
    exit(0);
  }
#endif

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

#ifdef DEBUG
  boil::oout<<"sharpen:iterate \n";
#endif

  /* iterate */
  for(int it=0; it<itmax; it++) {

#ifdef DEBUG
  boil::oout<<"it="<<it<<" "<<itmax<<"\n";
#endif

    stmp = 0.0;

    // x-direction
    //for_vmijk((*u),m,i,j,k){  //don't use vmijk. wall will be skipped!
    for(int i=si(); i<=ei()+1; i++)
      for(int j=sj(); j<=ej(); j++)
        for(int k=sk(); k<=ek(); k++){
#ifdef IB
          if (dom->ibody().off(Comp::u(),i,j,k)) continue;
#endif
          real tm, tc, fluxc, fluxd, coefc, coefd, ldt;
          // compression
          tm = sca[i-1][j][k] * (1.0-sca[i-1][j][k]) * nx[i-1][j][k];
          tc = sca[i]  [j][k] * (1.0-sca[i]  [j][k]) * nx[i]  [j][k];
          fluxc = 0.5 * (tm+tc) * dSx(i,j,k);
          coefc = beta(alp[i-1][j][k],alp[i][j][k],local)     // beta
                * 0.5*(sca.dxc(i-1)+sca.dxc(i))         // dx
                / (0.5*(atmp[i-1][j][k]+atmp[i][j][k]));// dxmax
          // diffusion
          fluxd = (sca[i][j][k]-sca[i-1][j][k])*dSx(i,j,k)/sca.dxw(i);
          coefd = beta(alp[i-1][j][k],alp[i][j][k],local)  // beta
                * epss * 0.5 * (sca.dxc(i-1)+sca.dxc(i));    // epss
          ldt  = 0.5*(fn [i-1][j][k]+fn [i][j][k]);

          // set ldt=0 instead of flux=0
          if(i==si()   && imin) ldt=0.0;
          if(i==ei()+1 && imax) ldt=0.0;
#ifdef WALL
          if(i==si()+1 && iminw) ldt=0.0;
          if(i==ei()   && imaxw) ldt=0.0;
#endif
#ifdef IB
          if(dom->ibody().off(Comp::u(),i-1,j,k)) ldt=0.0; // adjascent wall
          if(dom->ibody().off(Comp::u(),i+1,j,k)) ldt=0.0; // adjascent wall
#endif
      // add flux
      stmp[i-1][j][k] -= (fluxc*coefc - fluxd*coefd)*ldt;
      stmp[i  ][j][k] += (fluxc*coefc - fluxd*coefd)*ldt;
    }

#ifdef DEBUG
  boil::oout<<"y-direction, it="<<it<<" "<<itmax<<"\n";
#endif

    // y-direction
    //for_vmijk((*u),m,i,j,k){  //don't use vmijk. wall will be skipped!
    for(int i=si(); i<=ei(); i++)
      for(int j=sj(); j<=ej()+1; j++)
        for(int k=sk(); k<=ek(); k++){
#ifdef IB
          if (dom->ibody().off(Comp::v(),i,j,k)) continue;
#endif
          real tm, tc, fluxc, fluxd, coefc, coefd, ldt;
          // compression
          tm = sca[i][j-1][k] * (1.0-sca[i][j-1][k]) * ny[i][j-1][k];
          tc = sca[i][j]  [k] * (1.0-sca[i][j]  [k]) * ny[i][j]  [k];
          fluxc = 0.5 * (tm+tc) * dSy(i,j,k);
          coefc = beta(alp[i][j-1][k],alp[i][j][k],local)     // beta
                * 0.5*(sca.dyc(j-1)+sca.dyc(j))         // dy
                / (0.5*(atmp[i][j-1][k]+atmp[i][j][k]));// dxmax
          // diffusion
          fluxd = (sca[i][j][k]-sca[i][j-1][k])*dSy(i,j,k)/sca.dys(j);
          coefd = beta(alp[i][j-1][k],alp[i][j][k],local)  // beta
                * epss * 0.5 * (sca.dyc(j-1)+sca.dyc(j));    // epss
          ldt  = 0.5*(fn [i][j-1][k]+fn [i][j][k]);

          // set ldt=0 instead of flux=0
          if(j==sj()   && jmin) ldt=0.0;
          if(j==ej()+1 && jmax) ldt=0.0;
#ifdef WALL
          if(j==sj()+1 && jminw) ldt=0.0;
          if(j==ej()   && jmaxw) ldt=0.0;
#endif
#ifdef IB
          if(dom->ibody().off(Comp::v(),i,j-1,k)) ldt=0.0; // adjascent wall
          if(dom->ibody().off(Comp::v(),i,j+1,k)) ldt=0.0; // adjascent wall
#endif
      // add flux
      stmp[i][j-1][k] -= (fluxc*coefc - fluxd*coefd)*ldt;
      stmp[i][j  ][k] += (fluxc*coefc - fluxd*coefd)*ldt;
    }

    // z-direction
    //for_vmijk((*u),m,i,j,k){  //don't use vmijk. wall will be skipped!
    for(int i=si(); i<=ei(); i++)
      for(int j=sj(); j<=ej(); j++)
        for(int k=sk(); k<=ek()+1; k++){
#ifdef IB
          if (dom->ibody().off(Comp::w(),i,j,k)) continue;
#endif
          real tm, tc, fluxc, fluxd, coefc, coefd, ldt;
          // compression
          tm = sca[i][j][k-1] * (1.0-sca[i][j][k-1]) * nz[i][j][k-1];
          tc = sca[i][j][k]   * (1.0-sca[i][j][k]  ) * nz[i][j][k];
          fluxc = 0.5 * (tm+tc) * dSz(i,j,k);
          coefc = beta(alp[i][j][k-1],alp[i][j][k],local)     // beta
                * 0.5*(sca.dzc(k-1)+sca.dzc(k))         // dz
                / (0.5*(atmp[i][j][k-1]+atmp[i][j][k]));// dxmax
          // diffusion
          fluxd = (sca[i][j][k]-sca[i][j][k-1])*dSz(i,j,k)/sca.dzb(k);
          coefd = beta(alp[i][j][k-1],alp[i][j][k],local)  // beta
                * epss * 0.5 * (sca.dzc(k-1)+sca.dzc(k));    // epss
          ldt  = 0.5*(fn [i][j][k-1]+fn [i][j][k]);

          // set ldt=0 instead of flux=0
          if(k==sk()   && kmin) ldt=0.0;
          if(k==ek()+1 && kmax) ldt=0.0;
#ifdef WALL
          if(k==sk()+1 && kminw) ldt=0.0;
          if(k==ek()   && kmaxw) ldt=0.0;
#endif
#ifdef IB
          if(dom->ibody().off(Comp::w(),i,j,k-1)) ldt=0.0; // adjascent wall
          if(dom->ibody().off(Comp::w(),i,j,k+1)) ldt=0.0; // adjascent wall
#endif
          // add flux
          stmp[i][j][k-1] -= (fluxc*coefc - fluxd*coefd)*ldt;
          stmp[i][j][k  ] += (fluxc*coefc - fluxd*coefd)*ldt;
    }

    /* advance */
    dsca_max=0.0;  // L_infinity
    //real dsca_sum=0.0;  // L_infinity
    int imx,jmx,kmx;
    for_ijk(i,j,k){
      real dsca = 1.0/dV(i,j,k) * stmp[i][j][k];
      sca[i][j][k] += dsca;
      //dsca_sum += dsca;
      if(fabs(dsca)>dsca_max){
        dsca_max=fabs(dsca);
        imx=i;
        jmx=j;
        kmx=k;
      }
    }
    boil::cart.max_real(&dsca_max);

#ifdef DEBUG
  boil::oout<<"sharpen:boundary, it="<<it<<" "<<itmax<<"\n";
#endif

    /* boundary condition */
    //sca.bnd_update();
    bdcond(sca);
    sca.exchange_all();
#ifdef IB
  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
  ib_ext_scalar(sca);
  ib_bdcond(sca); 
#endif
  }

#if 1
  // for 1D 
  for_i(i){
    real phiave=0.0;
    int icell=0;
    for_jk(j,k){
      phiave+=sca[i][j][k];
      icell++;
    }
    phiave=phiave/real(icell);
    for_ajk(j,k){
      sca[i][j][k]=phiave;
    }
  }
#endif

  return true;
}

/******************************************************************************/
real CIPCSL2::beta(const real a1, const real a2, const bool local){
/***************************************************************************//**
* \brief  Calculate coefficient Beta using liner interpolation.
*            input: a1, a2
*            output: bout
*******************************************************************************/
  if (!local) {
    /* global sharpening */
    return 1.0;
  } else {
    /* local sharpening */
    real ave=0.5*(a1+a2);
    real bout;
    const real bmin=0.01;
    const real bmax=1.0;
    const real alpsat=2.0;
    real grd=(bmax-bmin)/alpsat;
    if(ave<=0.0){
      bout = bmin;
    } else {
      bout = std::min(1.0,(ave-alpsat)*grd+1.0);
    }
    return bout;
  }
}
