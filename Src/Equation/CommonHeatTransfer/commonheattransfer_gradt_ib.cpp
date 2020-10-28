#include "commonheattransfer.h"

/******************************************************************************/
real CommonHeatTransfer::gradt_ib(const Sign dir, const Comp & mcomp,
                                  const int i, const int j, const int k,
                                  const Old old,
                                  Scalar & val) const {
/***************************************************************************//**
*  \brief calculate gradient of temperature near immersed bodies
*******************************************************************************/

  /* if ibody can co-exist with no conduction in solid, this needs to change */
  if(!solid())
    return 0.;

  int of(+1);
  if(dir<0)
    of = -1;

  real len_f, len_s;
  real lam_f, lam_s;
  real tmp_f, tmp_s;

  tmp_f = val[i][j][k];

  /* assuming zero eddy viscosity in the boundary cell...*/
  lam_f = (old==Old::yes) ? lambda_old(i,j,k) : lambda(i,j,k);
  
  if       (mcomp==Comp::i()) {
    len_f = 0.5*val.dxc(i);
    len_s = 0.5*val.dxc(i+of);
    lam_s = solid()->lambda(i+of,j,k);
    tmp_s = val[i+of][j][k];
  } else if(mcomp==Comp::j()) {
    len_f = 0.5*val.dyc(j);
    len_s = 0.5*val.dyc(j+of);
    lam_s = solid()->lambda(i,j+of,k);
    tmp_s = val[i][j+of][k];
  } else {
    len_f = 0.5*val.dzc(k);
    len_s = 0.5*val.dzc(k+of);
    lam_s = solid()->lambda(i,j,k+of);
    tmp_s = val[i][j][k+of];
  }
 
  real tmp_node = temperature_node(
                      htwallmodel->dirac_wall_source,
                      len_s/lam_s, tmp_s, 
                      len_f/lam_f, tmp_f);

  if(dir<0) /* fluid is above */
    return (tmp_f-tmp_node)/len_f;
  else
    return (tmp_node-tmp_f)/len_f;
}
