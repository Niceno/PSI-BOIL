#include "phasechange.h"
//#define ST_LENGTH
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange::dmicro_intermediate(){

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=dom->ibody().nwx(i,j,k);
    real uy=dom->ibody().nwy(i,j,k);
    real uz=dom->ibody().nwz(i,j,k);
    Dir d = Dir::undefined();
    if (fabs(uz)>0.707) {
      d = Dir::kmin();
      if(uz>0.707){
        d = Dir::kmax();
      }
    } else if (fabs(ux)>0.707) {
      d = Dir::imin();
      if(ux>0.707){
        d = Dir::imax();
      }
    } else if (fabs(uy)>0.707) {
      d = Dir::jmin();
      if(uy>0.707){
        d = Dir::jmax();
      }
    } else {
      std::cout<<"phasechange_micro: Underdevelopment!!!\n";
      exit(0);
    }

    int iof=0, jof=0, kof=0;
    int ndir = 6;  // ibody()
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*--------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in fluid domain  |
    |         (i+iof, j+jof, k+kof) is in solid domain  |
    +--------------------------------------------------*/
    int ii=i+iof;
    int jj=j+jof;
    int kk=k+kof;

    /* include vapor */
    real a_vapor=nucl->area_vapor(i,j,k,d);
    if ( a_vapor > 0.0) {

      /* initialize */
      if ( nucl->dmicro[i][j][k] > boil::mega ) {
        nucl->dmicro[i][j][k] = nucl->dmicro0(i,j,k);
      }

      /* interface cell */
      if (nucl->dSprev[i][j] < clr.dSz(i,j,k)) {
        real s_prev = nucl->dSprev[i][j];
        real s_now  = a_vapor;
        if ( (s_now > s_prev) &&
             (nucl->dmicro[i][j][k] > nucl->dmicro_min+boil::pico) &&
             (nucl->dmicro[i][j][k] < boil::mega)) {
          real d_new = ( nucl->dmicro[i][j][k] * s_prev
                       + nucl->dmicro0(i,j,k) * (s_now - s_prev) )/s_now;
          nucl->dmicro[i][j][k] = max(d_new, nucl->dmicro_min);
        }
      }
    }
  }
#endif
}
