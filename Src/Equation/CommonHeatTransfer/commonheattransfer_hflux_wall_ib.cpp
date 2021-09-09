#include "commonheattransfer.h"

/***************************************************************************//**
*  \brief calculate heat flux on wall
*******************************************************************************/
real CommonHeatTransfer::hflux_wall_ib(const Scalar & val,
                                       const Scalar * diff_eddy) const {

  /* if ibody can co-exist with no conduction in solid, this needs to change */
  if(!solid())
    return 0.;

  real hflux=0.0;
  real areaw=0.0;

  for(int cc=0; cc<val.domain()->ibody().nccells(); cc++){
    int i,j,k;
    val.domain()->ibody().ijk(cc,&i,&j,&k);

    real ts, tw;  /* temperature of solid, wall */
    real ls;      /* lambda */
    real fd, area, gradt;

    /* w */
    if(val.domain()->ibody().off(i-1,j,k)){
      ts = val[i-1][j][k];
      tw = bndtpr_sol[Comp::i()][i][j][k];

      ls = solid()->lambda(i-1,j,k);
      fd = val.domain()->ibody().fdxw(i,j,k);
      real dxs = (1.0-fd)*val.dxw(i);
      area = val.dSx(Sign::neg(),i,j,k);
      real resists = dxs/ls;

      gradt = (ts-tw)/resists;
      hflux += area*(gradt+dirac_wall_source(i,j,k));
    }
    /* e */
    if(val.domain()->ibody().off(i+1,j,k)){
      ts = val[i+1][j][k];
      tw = bndtpr_sol[Comp::i()][i+1][j][k];

      ls = solid()->lambda(i+1,j,k);
      fd = val.domain()->ibody().fdxe(i,j,k);
      real dxs = (1.0-fd)*val.dxe(i);
      area = val.dSx(Sign::pos(),i,j,k);
      real resists = dxs/ls;

      gradt = (ts-tw)/resists;
      hflux += area*(gradt+dirac_wall_source(i,j,k));
    }
    /* s */
    if(val.domain()->ibody().off(i,j-1,k)){
      ts = val[i][j-1][k];
      tw = bndtpr_sol[Comp::j()][i][j][k];

      ls = solid()->lambda(i,j-1,k);
      fd = val.domain()->ibody().fdys(i,j,k);
      real dys = (1.0-fd)*val.dys(j);
      area = val.dSy(Sign::neg(),i,j,k);
      real resists = dys/ls;

      gradt = (ts-tw)/resists;
      hflux += area*(gradt+dirac_wall_source(i,j,k));
    }
    /* n */
    if(val.domain()->ibody().off(i,j+1,k)){
      ts = val[i][j+1][k];
      tw = bndtpr_sol[Comp::j()][i][j+1][k];

      ls = solid()->lambda(i,j+1,k);
      fd = val.domain()->ibody().fdyn(i,j,k);
      real dys = (1.0-fd)*val.dyn(j);
      area = val.dSy(Sign::pos(),i,j,k);
      real resists = dys/ls;

      gradt = (ts-tw)/resists;
      hflux += area*(gradt+dirac_wall_source(i,j,k));
    }
    /* b */
    if(val.domain()->ibody().off(i,j,k-1)){
      ts = val[i][j][k-1];
      tw = bndtpr_sol[Comp::k()][i][j][k];

      ls = solid()->lambda(i,j,k-1);
      fd = val.domain()->ibody().fdzb(i,j,k);
      real dzs = (1.0-fd)*val.dzb(k);
      area = val.dSz(Sign::neg(),i,j,k);
      real resists = dzs/ls;

      gradt = (ts-tw)/resists;
      hflux += area*(gradt+dirac_wall_source(i,j,k));
    }
    /* t */
    if(val.domain()->ibody().off(i,j,k+1)){
      ts = val[i][j][k+1];
      tw = bndtpr_sol[Comp::k()][i][j][k+1];

      ls = solid()->lambda(i,j,k+1);
      fd = val.domain()->ibody().fdzt(i,j,k);
      real dzs = (1.0-fd)*val.dzt(k);
      area = val.dSz(Sign::pos(),i,j,k);
      real resists = dzs/ls;

      gradt = (ts-tw)/resists;
      hflux += area*(gradt+dirac_wall_source(i,j,k));
    }
    areaw += area;
  }

  boil::cart.sum_real(&hflux);
  boil::cart.sum_real(&areaw);

  //boil::oout<<"areaw= "<<areaw<<"\n";

  return hflux;
}	

