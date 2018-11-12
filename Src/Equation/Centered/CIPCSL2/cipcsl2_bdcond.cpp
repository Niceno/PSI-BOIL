#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::bdcond(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for color function
******************************************************************************/

  Formula F;

  int i,j,k;

  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    // set loop range
    int ist,ied,jst,jed,kst,ked;
    if(sca.bc().at(b).si()==sca.bc().at(b).ei()){
      ist=ied=sca.bc().at(b).si();
    } else {
      if(sca.bc().at(b).si()==si()){ist=si()-1;}
      else{ist=sca.bc().at(b).si();}
      if(sca.bc().at(b).ei()==ei()){ied=ei()+1;}
      else{ied=sca.bc().at(b).ei();}
    }
    if(sca.bc().at(b).sj()==sca.bc().at(b).ej()){
      jst=jed=sca.bc().at(b).sj();
    } else {
      if(sca.bc().at(b).sj()==sj()){jst=sj()-1;}
      else{jst=sca.bc().at(b).sj();}
      if(sca.bc().at(b).ej()==ej()){jed=ej()+1;}
      else{jed=sca.bc().at(b).ej();}
    }
    if(sca.bc().at(b).sk()==sca.bc().at(b).ek()){
      kst=ked=sca.bc().at(b).sk();
    } else {
      if(sca.bc().at(b).sk()==sk()){kst=sk()-1;}
      else{kst=sca.bc().at(b).sk();}
      if(sca.bc().at(b).ek()==ek()){ked=ek()+1;}
      else{ked=sca.bc().at(b).ek();}
    }
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
      }
      /* formula is not defined */
      else {
        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(i=ist; i<=ied; i++)
          for(j=jst; j<=jed; j++)
            for(k=kst; k<=ked; k++) {
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
        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(i=ist; i<=ied; i++)
          for(j=jst; j<=jed; j++)
            for(k=kst; k<=ked; k++)
              sca[i][j][k]=std::max(0.0,min(1.0,sca[i+iof][j+jof][k+kof]));
      }
    }
  }
}
