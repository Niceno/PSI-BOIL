#include "vofaxisym.h"

/******************************************************************************/
real VOFaxisym::linf_scalar_error(const Scalar & sca, const Scalar & scb
                                  //,int & imax, int & jmax, int & kmax, bool & ebool
                                  ) {
/***************************************************************************//**
 \brief Calculate the Linf difference between two scalar fields.
    output: Linf
*******************************************************************************/

  real linferr(-1.0);
  for_vijk(sca,i,j,k) {
    if(dom->ibody().off(i,j,k)) continue;
    real valcc = sca[i  ][j][k  ];
    real valmc = sca[i-1][j][k  ];
    real valpc = sca[i+1][j][k  ];
    real valcm = sca[i  ][j][k-1];
    real valcp = sca[i  ][j][k+1];
    real valmm = sca[i-1][j][k-1];
    real valpm = sca[i+1][j][k-1];
    real valmp = sca[i-1][j][k+1];
    real valpp = sca[i+1][j][k+1];

    bool elvibool = (valcc-phisurf)*(valmc-phisurf)<=0.0
     ||(valcc-phisurf)*(valpc-phisurf)<=0.0
     ||(valcc-phisurf)*(valcm-phisurf)<=0.0
     ||(valcc-phisurf)*(valcp-phisurf)<=0.0
     ||(valcc-phisurf)*(valmm-phisurf)<=0.0
     ||(valcc-phisurf)*(valpm-phisurf)<=0.0
     ||(valcc-phisurf)*(valmp-phisurf)<=0.0
     ||(valcc-phisurf)*(valpp-phisurf)<=0.0;

    if(!elvibool) continue;

    real sdiff = fabs(std::max(0.0,std::min(1.0,sca[i][j][k]))
                     -std::max(0.0,std::min(1.0,scb[i][j][k])));
    if(sdiff>linferr) {
      linferr = sdiff;
#if 0
      imax = i;
      jmax = j;
      kmax = k;
      ebool = elvibool;
#endif
    }
  }
  boil::cart.max_real(&linferr);

  return linferr;
}
