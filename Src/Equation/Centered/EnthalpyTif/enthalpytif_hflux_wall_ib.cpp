#include "enthalpytif.h"
using namespace std;

/***************************************************************************//**
*  \brief calculate heat flux on wall
*******************************************************************************/
real EnthalpyTIF::hflux_wall_ib(const Scalar * diff_eddy) {

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
    if ((*clr)[i][j][k]<clrsurf) {
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
      tf = phi[i][j][k];
      ts = phi[i-1][j][k];
      ls = solid()->lambda(i-1,j,k);
      fd = dom->ibody().fdxw(i,j,k);
      real dxs = (1.0-fd)*phi.dxw(i);
      real dxf = fd*phi.dxw(i);
      area = phi.dSx(i,j,k);
      if(fs&&Interface(-1,Comp::i(),i,j,k)) {
        dxf = 0.5*phi.dxc(i) - distance_x(i,j,k,-1,tf);
        /* inversion of lambda */
        if((*clr)[i][j][k]<clrsurf) {
          lf = lambdal;
          cp_mass = cpl/rhol;
        } else {
          lf = lambdav;
          cp_mass = cpv/rhov;
        }
        if(diff_eddy){
          lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
        }
      }
      tw = temperature_node(dxs, ls, ts, dxf, lf, tf);
      gradt = (ts-tw)/dxs;
      hflux += area*ls*gradt;
    }
    // e
    if(dom->ibody().off(i+1,j,k)){
      tf = phi[i][j][k];
      ts = phi[i+1][j][k];
      ls = solid()->lambda(i+1,j,k);
      fd = dom->ibody().fdxe(i,j,k);
      real dxs = (1.0-fd)*phi.dxe(i);
      real dxf = fd*phi.dxe(i);
      area = phi.dSx(i,j,k);
      if(fs&&Interface(+1,Comp::i(),i,j,k)) {
        dxf = 0.5*phi.dxc(i) - distance_x(i,j,k,+1,tf);
        /* inversion of lambda */
        if((*clr)[i][j][k]<clrsurf) {
          lf = lambdal;
          cp_mass = cpl/rhol;
        } else {
          lf = lambdav;
          cp_mass = cpv/rhov;
        }
        if(diff_eddy){
          lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
        }
      }
      tw = temperature_node(dxs, ls, ts, dxf, lf, tf);
      gradt = (ts-tw)/dxs;
      hflux += area*ls*gradt;
    }
    // s
    if(dom->ibody().off(i,j-1,k)){
      tf = phi[i][j][k];
      ts = phi[i][j-1][k];
      ls = solid()->lambda(i,j-1,k);
      fd = dom->ibody().fdys(i,j,k);
      real dys = (1.0-fd)*phi.dys(j);
      real dyf = fd*phi.dys(j);
      area = phi.dSy(i,j,k);
      if(fs&&Interface(-1,Comp::j(),i,j,k)) {
        dyf = 0.5*phi.dyc(j) - distance_y(i,j,k,-1,tf);
        /* inversion of lambda */
        if((*clr)[i][j][k]<clrsurf) {
          lf = lambdal;
          cp_mass = cpl/rhol;
        } else {
          lf = lambdav;
          cp_mass = cpv/rhov;
        }
        if(diff_eddy){
          lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
        }
      }
      tw = temperature_node(dys, ls, ts, dyf, lf, tf);
      gradt = (ts-tw)/dys;
      hflux += area*ls*gradt;
    }
    // n
    if(dom->ibody().off(i,j+1,k)){
      tf = phi[i][j][k];
      ts = phi[i][j+1][k];
      ls = solid()->lambda(i,j+1,k);
      fd = dom->ibody().fdyn(i,j,k);
      real dys = (1.0-fd)*phi.dyn(j);
      real dyf = fd*phi.dyn(j);
      area = phi.dSy(i,j,k);
      if(fs&&Interface(+1,Comp::j(),i,j,k)) {
        dyf = 0.5*phi.dyc(j) - distance_y(i,j,k,+1,tf);
        /* inversion of lambda */
        if((*clr)[i][j][k]<clrsurf) {
          lf = lambdal;
          cp_mass = cpl/rhol;
        } else {
          lf = lambdav;
          cp_mass = cpv/rhov;
        }
        if(diff_eddy){
          lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
        }
      }
      tw = temperature_node(dys, ls, ts, dyf, lf, tf);
      gradt = (ts-tw)/dys;
      hflux += area*ls*gradt;
    }
    /* b */
    if(dom->ibody().off(i,j,k-1)){
      tf = phi[i][j][k];
      ts = phi[i][j][k-1];
      ls = solid()->lambda(i,j,k-1);
      fd = dom->ibody().fdzb(i,j,k);
      real dzs = (1.0-fd)*phi.dzb(k);
      real dzf = fd*phi.dzb(k);
      area = phi.dSz(i,j,k);
      if(fs&&Interface(-1,Comp::k(),i,j,k)) {
        dzf = 0.5*phi.dzc(k) - distance_z(i,j,k,-1,tf);
        /* inversion of lambda */
        if((*clr)[i][j][k]<clrsurf) {
          lf = lambdal;
          cp_mass = cpl/rhol;
        } else {
          lf = lambdav;
          cp_mass = cpv/rhov;
        }
        if(diff_eddy){
          lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
        }
      }
      tw = temperature_node(dzs, ls, ts, dzf, lf, tf);
      gradt = (ts-tw)/dzs;
      hflux += area*ls*gradt;
    }
    /* t */
    if(dom->ibody().off(i,j,k+1)){
      tf = phi[i][j][k];
      ts = phi[i][j][k+1];
      ls = solid()->lambda(i,j,k+1);
      fd = dom->ibody().fdzt(i,j,k);
      real dzs = (1.0-fd)*phi.dzt(k);
      real dzf = fd*phi.dzt(k);
      area = phi.dSz(i,j,k);
      if(fs&&Interface(+1,Comp::k(),i,j,k)) {
        dzf = 0.5*phi.dzc(k) - distance_z(i,j,k,+1,tf);
        /* inversion of lambda */
        if((*clr)[i][j][k]<clrsurf) {
          lf = lambdal;
          cp_mass = cpl/rhol;
        } else {
          lf = lambdav;
          cp_mass = cpv/rhov;
        }
        if(diff_eddy){
          lf += (*diff_eddy)[i][j][k]*cp_mass/turbP;
        }
      }
      tw = temperature_node(dzs, ls, ts, dzf, lf, tf);
      gradt = (ts-tw)/dzs;
      hflux += area*ls*gradt;
    }
    areaw += area;
  }

  boil::cart.sum_real(&hflux);
  boil::cart.sum_real(&areaw);

  //boil::oout<<"areaw= "<<areaw<<"\n";

  return hflux;
}	

