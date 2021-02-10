#include "commonheattransfer.h"

/******************************************************************************/
real CommonHeatTransfer::gradt_ib(const Sign dir, const Comp & mcomp,
                                  const int i, const int j, const int k,
                                  const Old old,
                                  Scalar & val) const {
/***************************************************************************//**
*  \brief calculate gradient of temperature near immersed bodies
*         - used by convection-old
*******************************************************************************/

  /* if ibody can co-exist with no conduction in solid, this needs to change */
  if(!solid())
    return 0.;

  int off(+1);
  //int of(+1), off(+1);
  if(dir<0) {
    //of = -1;
    off = 0;
  }

  real tmp_f = val[i][j][k];
  real tmp_node = bndtpr_flu[mcomp][i+off*(mcomp==Comp::i())]
                                   [j+off*(mcomp==Comp::j())]
                                   [k+off*(mcomp==Comp::k())];

  real len_f = distance_face(dir,mcomp,i,j,k);

  if(dir<0) /* fluid is above */
    return (tmp_f-tmp_node)/len_f;
  else
    return (tmp_node-tmp_f)/len_f;
}
