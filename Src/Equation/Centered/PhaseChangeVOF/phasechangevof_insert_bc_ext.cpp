#include "phasechangevof.h"
using namespace std;

/* obsolete!!!! */
/******************************************************************************/
void PhaseChangeVOF::insert_bc_ext(const Comp mcomp) {
/***************************************************************************//**
*  \brief correct extrapolation flag at walls
*  Cells with subgrid interfaces will not be extrapolated, since the gradt vals
*  - are computed for the fluid in which the cell centre is
*  - are computed for the fluid which is subgrid
*  I.e. both values are computed
*******************************************************************************/
  Dir dm = Dir::undefined();
  Dir dp = Dir::undefined();
  int dirm(-1), ofxm(0), ofym(0), ofzm(0);
  int dirp(+1), ofxp(0), ofyp(0), ofzp(0);

  if       (mcomp==Comp::i()) {
    dm = Dir::imin();
    dp = Dir::imax();
    ofxp = -1;
    ofxm = +1;
  } else if(mcomp==Comp::j()) {
    dm = Dir::jmin();
    dp = Dir::jmax();
    ofyp = -1;
    ofym = +1;
  } else {
    dm = Dir::kmin();
    dp = Dir::kmax();
    ofzp = -1;
    ofzm = +1;
  }

  for(int b = 0; b < clr.bc().count(); b++) {

    if(clr.bc().type_decomp(b))
      continue;

    if(clr.bc().type(b) == BndType::wall()) {

      Dir d = clr.bc().direction(b);
      if       (d == dm) {
        for_vijk( clr.bc().at(b), i,j,k ) {
          int ii = i+ofxm;
          int jj = j+ofym;
          int kk = k+ofzm;
          /* is there an interface between cell centre and wall? */
          if(Interface(dirm,mcomp,ii,jj,kk)) {
            stmp[ii][jj][kk] = 1; /* will not be extrapolated */
          }
        }
      } else if(d == dp) {
        for_vijk( clr.bc().at(b), i,j,k ) {
          int ii = i+ofxp;
          int jj = j+ofyp;
          int kk = k+ofzp;
          /* is there an interface between cell centre and wall? */
          if(Interface(dirp,mcomp,ii,jj,kk)) {
            stmp[ii][jj][kk] = 1; /* will not be extrapolated */
          }
        }
      } 

    } /* is wall? */
  } /* loop over bcs */

  return;
}
