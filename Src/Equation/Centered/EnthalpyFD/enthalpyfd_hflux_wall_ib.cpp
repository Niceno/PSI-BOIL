#include "enthalpyfd.h"
using namespace std;

/***************************************************************************//**
*  \brief calculate heat flux on immersed boundary in entire domain
*******************************************************************************/
void EnthalpyFD::hflux_wall_ib(const Scalar * diff_eddy) {
  hflux_wall_ib(Range<real>(-boil::yotta, boil::yotta),
                Range<real>(-boil::yotta, boil::yotta), diff_eddy);
  return;
}

/***************************************************************************//**
*  \brief calculate heat flux on immersed boundary
*  in restricted range
*******************************************************************************/
void EnthalpyFD::hflux_wall_ib(Range<real> xr, Range<real> yr,
                               const Scalar * diff_eddy) {

  // initialize
  real hflux[2] = {0.0,0.0};
  real area_sum[2] = {0.0,0.0};

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    real tf, ts, tw;  // temperature of fluid, solid, wall
    real lf, ls;      // lambda
    real fd, area, gradt;

    tf = phi[i][j][k];
    real cp_mass;
    if ((*clr)[i][j][k]<0.5) {
      lf = lambdav;
      cp_mass = cpv/rhov;
    } else {
      lf = lambdal;
      cp_mass = cpl/rhol;
    }
    if(diff_eddy){
      lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
    }

    // w
    if(dom->ibody().off(i-1,j,k)){
      ts = phi[i-1][j][k];
      ls = solid()->lambda(i-1,j,k);
      fd = dom->ibody().fdxw(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSx(i,j,k);
      gradt = (ts-tw)/(phi.dxw(i)*(1.0-fd));
    }
    // e
    if(dom->ibody().off(i+1,j,k)){
      ts = phi[i+1][j][k];
      ls = solid()->lambda(i+1,j,k);
      fd = dom->ibody().fdxe(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSx(i,j,k);
      gradt = (ts-tw)/(phi.dxe(i)*(1.0-fd));
    }
    // s
    if(dom->ibody().off(i,j-1,k)){
      ts = phi[i][j-1][k];
      ls = solid()->lambda(i,j-1,k);
      fd = dom->ibody().fdys(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSy(i,j,k);
      gradt = (ts-tw)/(phi.dys(j)*(1.0-fd));
    }
    // n
    if(dom->ibody().off(i,j+1,k)){
      ts = phi[i][j+1][k];
      ls = solid()->lambda(i,j+1,k);
      fd = dom->ibody().fdyn(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSy(i,j,k);
      gradt = (ts-tw)/(phi.dyn(j)*(1.0-fd));
    }
    /* b */
    if(dom->ibody().off(i,j,k-1)){
      ts = phi[i][j][k-1];
      ls = solid()->lambda(i,j,k-1);
      fd = dom->ibody().fdzb(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSz(i,j,k);
      gradt = (ts-tw)/(phi.dzb(k)*(1.0-fd));
    }
    /* t */
    if(dom->ibody().off(i,j,k+1)){
      ts = phi[i][j][k+1];
      ls = solid()->lambda(i,j,k+1);
      fd = dom->ibody().fdzt(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSz(i,j,k);
      gradt = (ts-tw)/(phi.dzt(k)*(1.0-fd));
    }

    if (xr.contains((*clr).xc(i)) && yr.contains((*clr).yc(j)) ) { // crude code
      if ((*clr)[i][j][k]>0.5) {
        hflux[0] += area*ls*gradt;
        area_sum[0] += area;
      } else {
        hflux[1] += area*ls*gradt;
        area_sum[1] += area;
      }
    }
  }

  boil::cart.sum_real_n(&hflux[0],2);
  boil::cart.sum_real_n(&area_sum[0],2);

  boil::oout<<"enthalpyFD:hflux= "<<time->current_time()
            <<" hflux_liquid[W] "<<hflux[0]
            <<" hflux_vapor[W] "<<hflux[1]
            <<" area_liquid[m2] "<<area_sum[0]
            <<" area_vapord[m2] "<<area_sum[1]<<"\n";

  return;
}
