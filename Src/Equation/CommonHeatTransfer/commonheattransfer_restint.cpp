#include "commonheattransfer.h"

/***************************************************************************//** 
*  \brief calculate effect of interfacial heat transfer resistance
*******************************************************************************/
/* 
 * only liquid is considered for the moment
 *
 *
 * set of coords: cell with the interface, i.e. to
 * if cell_marker is > 0, material properties are of the other phase

table:

  .--------.-------.-------.
  | INT\CM |   +   |   -   |
  .--------.-------.-------.
  |   +    |   0   |   1   |
  .--------.-------.-------.
  |   -    |   1   |   0   |
  .--------.-------.-------.

  This convoluted way is used because I wanted to avoid checking if
  a cell, which is solid or in wall, is liquid/gas.
  Although update_at_walls should allow such a check, it is perhaps
  dangerous to rely on that. With the approach coded below, status of
  fluid cell is always checked.

*/

void CommonHeatTransfer::resTint(const Sign & dir, const Comp & m,
                                 const int i0, const int j0, const int k0,
                                 const int ii, const int ji, const int ki,
                                 const int i1, const int j1, const int k1,
                                 const real dist, const Sign & cell_marker,
                                 real & tint, const Old old) const {

  /* sanity check */
  assert(topo->domain()->ibody().on(i0,j0,k0));
  
  Sign phase_marker = cell_marker*topo->sign_interface(ii,ji,ki,old);
  /* see explanation above */
  if(phase_marker < 0) {

    /* prepare stencil for evaluation */
    std::vector<StencilPoint> stencil;
    std::vector<real> coefs;
    stencil.push_back(StencilPoint(0,tint,0.));
    stencil.push_back(StencilPoint(1,tpr[i0][j0][k0],dist));

    /* resistance ghost distance */
    real res = int_resistance_liq(ii,ji,ki)*lambdal(ii,ji,ki);
    real resinv;
    if(m==Comp::i()) {
      resinv = fabs(topo->get_nx()[ii][ji][ki])/res;
    } else if(m==Comp::j()) {
      resinv = fabs(topo->get_ny()[ii][ji][ki])/res;
    } else {
      resinv = fabs(topo->get_nz()[ii][ji][ki])/res;
    }

    /* is there interface to the other side? */
    if(interface(dir,m,i0,j0,k0,old)) {

      /* default to first order for simplicity */
      topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::First());
      tint = (coefs[1]*stencil[1].val + resinv*tint)/(resinv - coefs[0]);
      return;

    /* are we at a solid-fluid boundary? */
    } else if(topo->domain()->ibody().off(i1,j1,k1)) {
      real dist1 = dist + distance_face(dir,m,i0,j0,k0);
      /* directional choice */
      if(dir>0) {
        stencil.push_back(StencilPoint(2,bndtpr_flu[m][i1][j1][k1],dist1));
      } else {
        stencil.push_back(StencilPoint(2,bndtpr_flu[m][i0][j0][k0],dist1));
      }

    /* neither interface, nor boundary */
    } else {
      real dist1 = dist + distance_center(dir,m,i0,j0,k0);
      stencil.push_back(StencilPoint(2,tpr[i1][j1][k1],dist1));
    }

    /* evaluate coefs */
    topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::Second());

    /* equation for resistance effect */
    tint = (coefs[1]*stencil[1].val + coefs[2]*stencil[2].val + resinv*tint)
         / (resinv - coefs[0]);

  } /* are we in liquid? */

  return;
}
