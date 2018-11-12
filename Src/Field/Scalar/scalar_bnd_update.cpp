#include "scalar.h"

/******************************************************************************/
void Scalar::bnd_update() {
/*-----------------------------------------------------------------------------+
|  updates boundary condition values for a scalar variable.                    |
|                                                                              |
|  it doesn't work for convective boundary conditions, because the necessary   |
|  information on physical properties and transfer coeffients are missing.     |
|                                                                              |
|  scalar_exchange_all should take account of periodic condition               |
|    1st: phi.bnd_update();                                                    |
|    2nd: phi.exchange_all();                                                  |
+-----------------------------------------------------------------------------*/

  Formula F;

  int i,j,k;

  for( int b=0; b<bc().count(); b++ ) {

    if( bc().type_decomp(b) ) continue;

    // set loop range
    int ist,ied,jst,jed,kst,ked;
    if(bc().at(b).si()==bc().at(b).ei()){
      ist=ied=bc().at(b).si();
    } else {
      if(bc().at(b).si()==si()){ist=si()-1;}
      else{ist=bc().at(b).si();}
      if(bc().at(b).ei()==ei()){ied=ei()+1;}
      else{ied=bc().at(b).ei();}
    }
    if(bc().at(b).sj()==bc().at(b).ej()){
      jst=jed=bc().at(b).sj();
    } else {
      if(bc().at(b).sj()==sj()){jst=sj()-1;}
      else{jst=bc().at(b).sj();}
      if(bc().at(b).ej()==ej()){jed=ej()+1;}
      else{jed=bc().at(b).ej();}
    }
    if(bc().at(b).sk()==bc().at(b).ek()){
      kst=ked=bc().at(b).sk();
    } else {
      if(bc().at(b).sk()==sk()){kst=sk()-1;}
      else{kst=bc().at(b).sk();}
      if(bc().at(b).ek()==ek()){ked=ek()+1;}
      else{ked=bc().at(b).ek();}
    }

    /*========================+ 
    |  dirichlet (and inlet)  |
    +========================*/
    if( bc().type(b) == BndType::dirichlet() ||
        bc().type(b) == BndType::inlet() ) {

      /* formula is defined */
      if( bc().formula(b) ) {
        for_vijk( bc().at(b), i,j,k ) {
          std::stringstream x, y, z, f;
          x << "x=" << xc(i); F.evaluate(x);
          y << "y=" << yc(j); F.evaluate(y);
          z << "z=" << zc(k); F.evaluate(z);
          f << bc().formula(b);

          val[i][j][k] = F.evaluate(f);
        }
      }
      /* formula is not defined */
      else {
        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(i=ist; i<=ied; i++)
          for(j=jst; j<=jed; j++)
            for(k=kst; k<=ked; k++) {
              val[i][j][k] = bc().value(b);
            } 
      }
    }

    /*==========+ 
    |  others   |
    +==========*/
    if( bc().type(b) == BndType::neumann()  ||
        bc().type(b) == BndType::symmetry() ||
        bc().type(b) == BndType::wall()     ||
        bc().type(b) == BndType::outlet() ) {

      int iof=0, jof=0, kof=0;

      Dir d = bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      
      if( d == Dir::ibody() ) {
        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) ) 
           val[i-1][j][k] = val[i][j][k];
          if( dom->ibody().on(i-1,j,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i-1][j][k];
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) ) 
           val[i+1][j][k] = val[i][j][k];
          if( dom->ibody().on(i+1,j,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i+1][j][k];
  

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) ) 
           val[i][j-1][k] = val[i][j][k];
          if( dom->ibody().on(i,j-1,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j-1][k];
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) ) 
           val[i][j+1][k] = val[i][j][k];
          if( dom->ibody().on(i,j+1,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j+1][k];
          
  
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) ) 
           val[i][j][k-1] = val[i][j][k];
          if( dom->ibody().on(i,j,k-1) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j][k-1];
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) ) 
           val[i][j][k+1] = val[i][j][k];
          if( dom->ibody().on(i,j,k+1) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j][k+1];
        }

      } else if(d != Dir::undefined()) {
        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(i=ist; i<=ied; i++)
          for(j=jst; j<=jed; j++)
            for(k=kst; k<=ked; k++)
              val[i][j][k]=val[i+iof][j+jof][k+kof];
      }
    }

    /*======================+ 
    |  convective boundary  |
    +======================*/
    if( bc().type(b) == BndType::convective() ) {
      for_vijk( bc().at(b), i,j,k ) {
        /* do nothing here, it is handled in the Centered class */
      }
    }

  }
}
