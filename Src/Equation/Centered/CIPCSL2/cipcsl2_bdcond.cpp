#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::bdcond(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for color function
******************************************************************************/

  Formula F;

  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    /*========================+ 
    |  dirichlet (and inlet)  |
    +========================*/
    if( sca.bc().type(b) == BndType::dirichlet() ||
        sca.bc().type(b) == BndType::inlet() ) {

      /* formula is defined */
      if( sca.bc().formula(b) ) {
        for_vijk( sca.bc().at(b), i,j,k ) {
          std::stringstream x, y, z, f;
          x << "x=" << sca.xc(i); F.evaluate(x);
          y << "y=" << sca.yc(j); F.evaluate(y);
          z << "z=" << sca.zc(k); F.evaluate(z);
          f << sca.bc().formula(b);

          sca[i][j][k] = F.evaluate(f);
        }
      } else {
        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for_vijk( sca.bc().at(b), i,j,k ){
          if (i<=si()-2) continue;
          if (i>=ei()+2) continue;
          if (j<=sj()-2) continue;
          if (j>=ej()+2) continue;
          if (k<=sk()-2) continue;
          if (k>=ek()+2) continue;
          sca[i][j][k] = sca.bc().value(b);
        }
      }
    }

    /*==========+ 
    |  others   |
    +==========*/
    if( sca.bc().type(b) == BndType::neumann()  ||
        sca.bc().type(b) == BndType::symmetry() ||
        sca.bc().type(b) == BndType::wall()     ||
        sca.bc().type(b) == BndType::outlet() ) {

      int iof=0, jof=0, kof=0;

      Dir d = sca.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
#if 0
      if( d == Dir::ibody() ) {
        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) )
           sca[i-1][j][k] = sca[i][j][k];
          if( dom->ibody().on(i-1,j,k) && dom->ibody().off(i,j,k) )
           sca[i][j][k] = sca[i-1][j][k];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) )
           sca[i+1][j][k] = sca[i][j][k];
          if( dom->ibody().on(i+1,j,k) && dom->ibody().off(i,j,k) )
           sca[i][j][k] = sca[i+1][j][k];


          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) )
           sca[i][j-1][k] = sca[i][j][k];
          if( dom->ibody().on(i,j-1,k) && dom->ibody().off(i,j,k) )
           sca[i][j][k] = sca[i][j-1][k];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) )
           sca[i][j+1][k] = sca[i][j][k];
          if( dom->ibody().on(i,j+1,k) && dom->ibody().off(i,j,k) )
           sca[i][j][k] = sca[i][j+1][k];


          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) )
           sca[i][j][k-1] = sca[i][j][k];
          if( dom->ibody().on(i,j,k-1) && dom->ibody().off(i,j,k) )
           sca[i][j][k] = sca[i][j][k-1];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) )
           sca[i][j][k+1] = sca[i][j][k];
          if( dom->ibody().on(i,j,k+1) && dom->ibody().off(i,j,k) )
           sca[i][j][k] = sca[i][j][k+1];
        }
      } else if(d != Dir::undefined()) {
#endif
      if(d != Dir::undefined()) {
        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        for_vijk( sca.bc().at(b), i,j,k ){
          if (i<=si()-2) continue;
          if (i>=ei()+2) continue;
          if (j<=sj()-2) continue;
          if (j>=ej()+2) continue;
          if (k<=sk()-2) continue;
          if (k>=ek()+2) continue;
              sca[i][j][k]=std::max(0.0,min(1.0,sca[i+iof][j+jof][k+kof]));
        }
      }
    }
  }
}
