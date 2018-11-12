#include "enthalpyfd.h"
using namespace std;

/***************************************************************************//**
*  \brief calculate heat flux on wall
*******************************************************************************/
real EnthalpyFD::hflux_wall_ib(const Scalar * diff_eddy) {

  real hflux=0.0;
  real areaw=0.0;

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
      hflux += area*ls*gradt;
    }
    // e
    if(dom->ibody().off(i+1,j,k)){
      ts = phi[i+1][j][k];
      ls = solid()->lambda(i+1,j,k);
      fd = dom->ibody().fdxe(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSx(i,j,k);
      gradt = (ts-tw)/(phi.dxe(i)*(1.0-fd));
      hflux += area*ls*gradt;
    }
    // s
    if(dom->ibody().off(i,j-1,k)){
      ts = phi[i][j-1][k];
      ls = solid()->lambda(i,j-1,k);
      fd = dom->ibody().fdys(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSy(i,j,k);
      gradt = (ts-tw)/(phi.dys(j)*(1.0-fd));
      hflux += area*ls*gradt;
    }
    // n
    if(dom->ibody().off(i,j+1,k)){
      ts = phi[i][j+1][k];
      ls = solid()->lambda(i,j+1,k);
      fd = dom->ibody().fdyn(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSy(i,j,k);
      gradt = (ts-tw)/(phi.dyn(j)*(1.0-fd));
      hflux += area*ls*gradt;
    }
    /* b */
    if(dom->ibody().off(i,j,k-1)){
      ts = phi[i][j][k-1];
      ls = solid()->lambda(i,j,k-1);
      fd = dom->ibody().fdzb(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSz(i,j,k);
      gradt = (ts-tw)/(phi.dzb(k)*(1.0-fd));
      hflux += area*ls*gradt;
    }
    /* t */
    if(dom->ibody().off(i,j,k+1)){
      ts = phi[i][j][k+1];
      ls = solid()->lambda(i,j,k+1);
      fd = dom->ibody().fdzt(i,j,k);
      tw = (ls*ts*fd + lf*tf*(1.0-fd))/(ls*fd + lf*(1.0-fd));
      area = phi.dSz(i,j,k);
      gradt = (ts-tw)/(phi.dzt(k)*(1.0-fd));
      hflux += area*ls*gradt;
    }
    areaw += area;
  }

  boil::cart.sum_real(&hflux);
  boil::cart.sum_real(&areaw);

  //boil::oout<<"areaw= "<<areaw<<"\n";

  return hflux;
}
