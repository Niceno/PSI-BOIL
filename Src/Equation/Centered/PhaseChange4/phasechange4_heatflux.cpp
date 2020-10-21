#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::heat_flux(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief Calculate heat flux in cells.
*******************************************************************************/

  /* these are heat flux */
  txl = tyl = tzl = 0.;
  txv = tyv = tzv = 0.;
  tnl = tnv = 0.;

  /* calculate temperature at solid-fluid boundaries */
  if(solid())
    calculate_node_temperature(diff_eddy);

  /* calculate heat flux */
  cal_hf(diff_eddy);
  dirac_source_terms(); /* source terms *in* the solid-fluid wall */

  txl.bnd_update();
  tyl.bnd_update();
  tzl.bnd_update();
  txv.bnd_update();
  tyv.bnd_update();
  tzv.bnd_update();

  /* insert special boundary conditions */
  insert_bc_hf(diff_eddy); /* dirichlet */

  txl.exchange_all();
  tyl.exchange_all();
  tzl.exchange_all();
  txv.exchange_all();
  tyv.exchange_all();
  tzv.exchange_all();

#if 0
  boil::plot->plot(*(topo->clr),txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  boil::plot->plot(*(topo->clr),txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  exit(0);
#endif

#if 1
  /* extrapolate heat flux, values in brackets indicate
     iflag values for extrapolated cells */

  std::set<int> pos_ext = {1,2};
  std::set<int> neg_ext = {-1,-2};

  if(use_unconditional_extrapolation) {
    pos_ext = {-1,1,2};
    neg_ext = {1,-1,-2};
  }

  topo->extrapolate(txv,Sign::pos(),pos_ext);
  topo->extrapolate(tyv,Sign::pos(),pos_ext);
  topo->extrapolate(tzv,Sign::pos(),pos_ext);
  topo->extrapolate(txl,Sign::neg(),neg_ext);
  topo->extrapolate(tyl,Sign::neg(),neg_ext);
  topo->extrapolate(tzl,Sign::neg(),neg_ext);
#endif

  /* calculate the normal component */
  for_aijk(i,j,k) {
    tnv[i][j][k] = txv[i][j][k]*nx[i][j][k]
                 + tyv[i][j][k]*ny[i][j][k]
                 + tzv[i][j][k]*nz[i][j][k];
    tnl[i][j][k] = txl[i][j][k]*nx[i][j][k]
                 + tyl[i][j][k]*ny[i][j][k]
                 + tzl[i][j][k]*nz[i][j][k];
  }

#if 0
  if(solid()) {
    /* z direction */
    Comp m = Comp::w();

    for_ijk(i,j,k) {
      /* bottom is in wall and this is an interfacial cell */
      if(dom->ibody().off(i,j,k-1) && interface(i,j,k)) {
        tnl[i][j][k] = -lambda(i,j,k-1)*(bndtpr[m][i][j][k]-tpr[i][j][k-1])/(0.5*phi.dzc(k));
        tnv[i][j][k] = 0.0;
      } /* 1 above solid */
    } /* ijk */
    tnl.exchange();
    tnv.exchange();
  }
#endif
  
#if 0
  boil::plot->plot(*(topo->clr),txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  boil::plot->plot(*(topo->clr),txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  boil::plot->plot(*(topo->clr),tnl,tnv,"clr-tnl-tnv",time->current_step());
  exit(0);
#endif

  return;
}

