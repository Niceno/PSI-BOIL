#include "commonheattransfer.h"

/***************************************************************************//**
*  \brief calculate heat flux on wall
*******************************************************************************/
real CommonHeatTransfer::hflux_wall_ib(const Scalar & val,
                                       const Scalar * diff_eddy) const {

  /* if ibody can co-exist with no conduction in solid, this needs to change */
  if(!solid())
    return 0;

  real hflux=0.0;
  real areaw=0.0;

  for(int cc=0; cc<val.domain()->ibody().nccells(); cc++){
    int i,j,k;
    val.domain()->ibody().ijk(cc,&i,&j,&k);

    real tf, ts, tw;  // temperature of fluid, solid, wall
    real lf, ls;      // lambda
    real fd, area, gradt;

    tf = val[i][j][k];
    lf = lambda(i,j,k,diff_eddy);

    // w
    if(val.domain()->ibody().off(i-1,j,k)){
      tf = val[i][j][k];
      ts = val[i-1][j][k];
      ls = solid()->lambda(i-1,j,k);
      fd = val.domain()->ibody().fdxw(i,j,k);
      real dxs = (1.0-fd)*val.dxw(i);
      real dxf = fd*val.dxw(i);
      area = val.dSx(Sign::neg(),i,j,k);
      real resistf = dxf/lf;
      real resists = dxs/ls;
      if(interface(Sign::neg(),Comp::i(),i,j,k)) {
        //dxf = 0.5*val.dxc(i) - distance_int_x(Sign::neg(),i,j,k,tf);
        dxf = distance_int_x(Sign::pos(),i-1,j,k,tf)-0.5*val.dxc(i-1);
        /* inversion of lambda */
        lf = lambda_inv(i,j,k,diff_eddy);
        resistf = dxf/lf + htwallmodel->near_wall_resist;
      }
      tw = temperature_node(htwallmodel->dirac_wall_source,
                            resists, ts, resistf, tf);
      gradt = (ts-tw)/dxs;
      hflux += area*ls*gradt;
    }
    // e
    if(val.domain()->ibody().off(i+1,j,k)){
      tf = val[i][j][k];
      ts = val[i+1][j][k];
      ls = solid()->lambda(i+1,j,k);
      fd = val.domain()->ibody().fdxe(i,j,k);
      real dxs = (1.0-fd)*val.dxe(i);
      real dxf = fd*val.dxe(i);
      area = val.dSx(Sign::pos(),i,j,k);
      real resistf = dxf/lf;
      real resists = dxs/ls;
      if(interface(Sign::pos(),Comp::i(),i,j,k)) {
        //dxf = 0.5*val.dxc(i) - distance_int_x(Sign::pos(),i,j,k,tf);
        dxf = distance_int_x(Sign::neg(),i+1,j,k,tf)-0.5*val.dxc(i+1);
        /* inversion of lambda */
        lf = lambda_inv(i,j,k,diff_eddy);
        resistf = dxf/lf + htwallmodel->near_wall_resist;
      }
      tw = temperature_node(htwallmodel->dirac_wall_source,
                            resists, ts, resistf, tf);
      gradt = (ts-tw)/dxs;
      hflux += area*ls*gradt;
    }
    // s
    if(val.domain()->ibody().off(i,j-1,k)){
      tf = val[i][j][k];
      ts = val[i][j-1][k];
      ls = solid()->lambda(i,j-1,k);
      fd = val.domain()->ibody().fdys(i,j,k);
      real dys = (1.0-fd)*val.dys(j);
      real dyf = fd*val.dys(j);
      area = val.dSy(Sign::neg(),i,j,k);
      real resistf = dyf/lf;
      real resists = dys/ls;
      if(interface(Sign::neg(),Comp::j(),i,j,k)) {
        //dyf = 0.5*val.dyc(j) - distance_int_y(Sign::neg(),i,j,k,tf);
        dyf = distance_int_y(Sign::pos(),i,j-1,k,tf)-0.5*val.dyc(j-1);
        /* inversion of lambda */
        lf = lambda_inv(i,j,k,diff_eddy);
        resistf = dyf/lf + htwallmodel->near_wall_resist;
      }
      tw = temperature_node(htwallmodel->dirac_wall_source,
                            resists, ts, resistf, tf);
      gradt = (ts-tw)/dys;
      hflux += area*ls*gradt;
    }
    // n
    if(val.domain()->ibody().off(i,j+1,k)){
      tf = val[i][j][k];
      ts = val[i][j+1][k];
      ls = solid()->lambda(i,j+1,k);
      fd = val.domain()->ibody().fdyn(i,j,k);
      real dys = (1.0-fd)*val.dyn(j);
      real dyf = fd*val.dyn(j);
      area = val.dSy(Sign::pos(),i,j,k);
      real resistf = dyf/lf;
      real resists = dys/ls;
      if(interface(Sign::pos(),Comp::j(),i,j,k)) {
        //dyf = 0.5*val.dyc(j) - distance_int_y(Sign::pos(),i,j,k,tf);
        dyf = distance_int_y(Sign::neg(),i,j+1,k,tf)-0.5*val.dyc(j+1);
        /* inversion of lambda */
        lf = lambda_inv(i,j,k,diff_eddy);
        resistf = dyf/lf + htwallmodel->near_wall_resist;
      }
      tw = temperature_node(htwallmodel->dirac_wall_source,
                            resists, ts, resistf, tf);
      gradt = (ts-tw)/dys;
      hflux += area*ls*gradt;
    }
    /* b */
    if(val.domain()->ibody().off(i,j,k-1)){
      tf = val[i][j][k];
      ts = val[i][j][k-1];
      ls = solid()->lambda(i,j,k-1);
      fd = val.domain()->ibody().fdzb(i,j,k);
      real dzs = (1.0-fd)*val.dzb(k);
      real dzf = fd*val.dzb(k);
      area = val.dSz(Sign::neg(),i,j,k);
      real resistf = dzf/lf;
      real resists = dzs/ls;
      if(interface(Sign::neg(),Comp::k(),i,j,k)) {
        //dzf = 0.5*val.dzc(k) - distance_int_z(Sign::neg(),i,j,k,tf);
        dzf = distance_int_z(Sign::pos(),i,j,k-1,tf)-0.5*val.dzc(k-1);
        /* inversion of lambda */
        lf = lambda_inv(i,j,k,diff_eddy);
        resistf = dzf/lf + htwallmodel->near_wall_resist;
      }
      tw = temperature_node(htwallmodel->dirac_wall_source,
                            resists, ts, resistf, tf);
      gradt = (ts-tw)/dzs;
      hflux += area*ls*gradt;
    }
    /* t */
    if(val.domain()->ibody().off(i,j,k+1)){
      tf = val[i][j][k];
      ts = val[i][j][k+1];
      ls = solid()->lambda(i,j,k+1);
      fd = val.domain()->ibody().fdzt(i,j,k);
      real dzs = (1.0-fd)*val.dzt(k);
      real dzf = fd*val.dzt(k);
      area = val.dSz(Sign::pos(),i,j,k);
      real resistf = dzf/lf;
      real resists = dzs/ls;
      if(interface(Sign::pos(),Comp::k(),i,j,k)) {
        //dzf = 0.5*val.dzc(k) - distance_int_z(Sign::pos(),i,j,k,tf);
        dzf = distance_int_z(Sign::neg(),i,j,k+1,tf)-0.5*val.dzc(k+1);
        /* inversion of lambda */
        lf = lambda_inv(i,j,k,diff_eddy);
        resistf = dzf/lf + htwallmodel->near_wall_resist;
      }
      tw = temperature_node(htwallmodel->dirac_wall_source,
                            resists, ts, resistf, tf);
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

