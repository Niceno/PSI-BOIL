#include "topology.h"

/******************************************************************************/
void Topology::cal_fs_interp(const Scalar & scp, Vector & fs,
                             const real tol_wall, const bool use_subgrid) {
/***************************************************************************//**
 \brief Calculate free-surface position between cell centers
    if there is no interface in the cell, unreal=yotta (=1e+24) is stored.
    output: fs
*******************************************************************************/

  /* fs is a result of calculations involving color
   * and cell dimensions. Thus, by taking care of properly exchanging
   * these variables, we can safely venture to calculate fs in buffer
   * cells without dealing with the exchange() method, which would 
   * not function properly for fs. So don't worry and let's loop over
   * (almost) all cells, including buffers except for the furthest one */
  int ibeg(scp.si()), iend(scp.ei()+1);
  int jbeg(scp.sj()), jend(scp.ej()+1);
  int kbeg(scp.sk()), kend(scp.ek()+1);

  if(bflag_struct.ifull) {
    if(bflag_struct.iminp) {
      ibeg -= boil::BW-1;
    }
    if(bflag_struct.imaxp) {
      iend += boil::BW-1;
    }
  }
  if(bflag_struct.jfull) {
    if(bflag_struct.jminp) {
      jbeg -= boil::BW-1;
    }
    if(bflag_struct.jmaxp) {
      jend += boil::BW-1;
    }
  }
  if(bflag_struct.kfull) {
    if(bflag_struct.kminp) {
      kbeg -= boil::BW-1;
    }
    if(bflag_struct.kmaxp) {
      kend += boil::BW-1;
    }
  }
        
  /* initialize */
  for_m(m)
    for_avmijk(fs,m,i,j,k)
      fs[m][i][j][k] = boil::unreal;

  Comp m;

  /******************************************
  *             x-direction                 *
  ******************************************/

  m = Comp::i();
  for(int i=ibeg; i<=iend; i++)
  for(int j=scp.sj(); j<=scp.ej(); j++)
  for(int k=scp.sk(); k<=scp.ek(); k++) {

    /* degenerate cases */
    if(domain()->ibody().off(i-1,j,k)&&domain()->ibody().off(i,j,k))
      continue;

    real clrw = scp[i-1][j][k];
    real clre = scp[i  ][j][k];

    if((clrw-clrsurf)*(clre-clrsurf)>0.0)
      continue;

    if(  (clrw<boil::pico&&clre-1.0>-boil::pico)
       ||(clrw-1.0>-boil::pico&&clre<boil::pico)) {
      fs[m][i][j][k] = scp.xn(i);
      continue;
    }
 
    if(fabs(clrw-clre)<boil::pico) {
      fs[m][i][j][k] = scp.xn(i);
      continue;
    }

    fs[m][i][j][k] =  scp.xc(i-1) + (clrsurf-clrw)/(clre-clrw)*scp.dxe(i-1);
  }

  /******************************************
  *             y-direction                 *
  ******************************************/

  m = Comp::j();
  for(int i=scp.si(); i<=scp.ei(); i++)
  for(int j=jbeg; j<=jend; j++)
  for(int k=scp.sk(); k<=scp.ek(); k++) {

    /* degenerate cases */
    if(domain()->ibody().off(i,j-1,k)&&domain()->ibody().off(i,j,k))
      continue;

    real clrs = scp[i][j-1][k];
    real clrn = scp[i][j  ][k];

    if((clrs-clrsurf)*(clrn-clrsurf)>0.0)
      continue;

    if(  (clrs<boil::pico&&clrn-1.0>-boil::pico)
       ||(clrs-1.0>-boil::pico&&clrn<boil::pico)) {
      fs[m][i][j][k] = scp.yn(j);
      continue;
    }

    if(fabs(clrs-clrn)<boil::pico) {
      fs[m][i][j][k] = scp.yn(j);
      continue;
    }

    fs[m][i][j][k] =  scp.yc(j-1) + (clrsurf-clrs)/(clrn-clrs)*scp.dyn(j-1);
  }

  /******************************************
  *             z-direction                 *
  ******************************************/

  m = Comp::k();
  for(int i=scp.si(); i<=scp.ei(); i++)
  for(int j=scp.sj(); j<=scp.ej(); j++)
  for(int k=kbeg; k<=kend; k++) {

    /* degenerate cases */
    if(domain()->ibody().off(i,j,k-1)&&domain()->ibody().off(i,j,k))
      continue;

    real clrb = scp[i][j][k-1];
    real clrt = scp[i][j][k  ];

    if((clrb-clrsurf)*(clrt-clrsurf)>0.0)
      continue;

    if(  (clrb<boil::pico&&clrt-1.0>-boil::pico)
       ||(clrb-1.0>-boil::pico&&clrt<boil::pico)) {
      fs[m][i][j][k] = scp.zn(k);
      continue;
    }

    if(fabs(clrb-clrt)<boil::pico) {
      fs[m][i][j][k] = scp.zn(k);
      continue;
    }

    fs[m][i][j][k] =  scp.zc(k-1) + (clrsurf-clrb)/(clrt-clrb)*scp.dzt(k-1);
  }

  /* correct at boundaries */
  if(use_subgrid)
    fs_bnd_geometric(scp,fs,tol_wall);
  else
    fs_bnd_nosubgrid(scp,fs,tol_wall);
  //fs.exchange_all();

  //boil::plot->plot(fs,scp, "fs-clr", 0);

  return;
}
