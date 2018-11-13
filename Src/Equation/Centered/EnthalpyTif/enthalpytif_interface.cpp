#include "enthalpytif.h"

/***************************************************************************//**
 *  Checks if the given cell is at an interface
******************************************************************************/

/******************************************\
| if combination color + indicator is used |
\******************************************/
bool EnthalpyTIF::Interface(const int i, const int j, const int k,
                            const Scalar & heaviside) {
  if(Interface(heaviside[i][j][k])) {
    return true;
  } 
  return false;
}

bool EnthalpyTIF::Interface(const real heavi) {
  return heavi > boil::atto && heavi < 1.0-boil::atto;
}

/******************************************\
|            if vof is used                |
\******************************************/
#if 0
bool EnthalpyTIF::Interface(const int i, const int j, const int k) {
  if(intflag[i][j][k]>0) {
    return true;
  }
  return false;
}

bool EnthalpyTIF::Interface_old(const int i, const int j, const int k) {
  if(intflagold[i][j][k]>0) {
    return true;
  }
  return false;
}

/* x-direction */
bool EnthalpyTIF::Interface1D_x(const int i, const int j, const int k) {

  real intm = (*fs)[Comp::i()][i  ][j][k];
  real intp = (*fs)[Comp::i()][i+1][j][k];
  /* is result inside of cell + in the correct direction? */
  real centrex = phi.xc(i);
  real edgep   = centrex+0.5*phi.dxc(i);
  real edgem   = centrex-0.5*phi.dxc(i);
  if(  (intm>=edgem&&intm<=centrex)
     ||(intp<=edgep&&intp>=centrex)) {
      return true;
  }

  return false;
}

/* y-direction */
bool EnthalpyTIF::Interface1D_y(const int i, const int j, const int k) {

  real intm = (*fs)[Comp::j()][i][j  ][k];
  real intp = (*fs)[Comp::j()][i][j+1][k];
  /* is result inside of cell + in the correct direction? */
  real centrey = phi.yc(j);
  real edgep   = centrey+0.5*phi.dyc(j);
  real edgem   = centrey-0.5*phi.dyc(j);
  if(  (intm>=edgem&&intm<=centrey)
     ||(intp<=edgep&&intp>=centrey)) {
      return true;
  }

  return false;
}

/* z-direction */
bool EnthalpyTIF::Interface1D_z(const int i, const int j, const int k) {

  real intm = (*fs)[Comp::k()][i][j][k  ];
  real intp = (*fs)[Comp::k()][i][j][k+1];
  /* is result inside of cell + in the correct direction? */
  real centrez = phi.zc(k);
  real edgep   = centrez+0.5*phi.dzc(k);
  real edgem   = centrez-0.5*phi.dzc(k);
  if(  (intm>=edgem&&intm<=centrey)
     ||(intp<=edgep&&intp>=centrey)) {
      return true;
  }

  return false;
}
#endif
