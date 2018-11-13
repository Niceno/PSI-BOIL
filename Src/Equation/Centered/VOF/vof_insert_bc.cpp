#include "vof.h"

/******************************************************************************/
void VOF::insert_bc(const Scalar & val) {
/***************************************************************************//**
*  \brief Set boundary value for a scalar variable.
*         Value at corner points, such as val[0][0][0], will be set as well.
*         scalar_exchange(_all) should take account of periodic condition.
*           1st: inset_bc(phi);
*           2nd: phi.exchange_all();
*******************************************************************************/

  Formula form;

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    /* set loop range */
    int ist,ied,jst,jed,kst,ked;
    if(val.bc().at(b).si()==val.bc().at(b).ei()){
      ist=ied=val.bc().at(b).si();
    } else {
      if(val.bc().at(b).si()==si()){ist=si()-1;}
      else{ist=val.bc().at(b).si();}
      if(val.bc().at(b).ei()==ei()){ied=ei()+1;}
      else{ied=val.bc().at(b).ei();}
    }
    if(val.bc().at(b).sj()==val.bc().at(b).ej()){
      jst=jed=val.bc().at(b).sj();
    } else {
      if(val.bc().at(b).sj()==sj()){jst=sj()-1;}
      else{jst=val.bc().at(b).sj();}
      if(val.bc().at(b).ej()==ej()){jed=ej()+1;}
      else{jed=val.bc().at(b).ej();}
    }
    if(val.bc().at(b).sk()==val.bc().at(b).ek()){
      kst=ked=val.bc().at(b).sk();
    } else {
      if(val.bc().at(b).sk()==sk()){kst=sk()-1;}
      else{kst=val.bc().at(b).sk();}
      if(val.bc().at(b).ek()==ek()){ked=ek()+1;}
      else{ked=val.bc().at(b).ek();}
    }

    /*------------------------+ 
    |  dirichlet (and inlet)  |
    +------------------------*/
    if( val.bc().type(b) == BndType::dirichlet() ||
        val.bc().type(b) == BndType::inlet() ) {

      /* formula is defined */
      if( val.bc().formula(b) ) {
        for_vijk( val.bc().at(b), i,j,k ) {
          std::stringstream x, y, z, f;
          x << "x=" << val.xc(i); form.evaluate(x);
          y << "y=" << val.yc(j); form.evaluate(y);
          z << "z=" << val.zc(k); form.evaluate(z);
          f << val.bc().formula(b);

          val[i][j][k] = form.evaluate(f);
        }
      }
      /* formula is not defined */
      else {
        //debug: for_vijk( val.bc().at(b), i,j,k )
        for(i=ist; i<=ied; i++)
        for(j=jst; j<=jed; j++)
        for(k=kst; k<=ked; k++)
          val[i][j][k] = val.bc().value(b);
      }
    }

    /*---------+ 
    |  others  |
    +---------*/
    if( val.bc().type(b) == BndType::neumann()
      ||val.bc().type(b) == BndType::symmetry() 
      ||val.bc().type(b) == BndType::wall()
      ||val.bc().type(b) == BndType::outlet() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {

        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        //debug: for_vijk( val.bc().at(b), i,j,k )
        for(i=ist; i<=ied; i++)
        for(j=jst; j<=jed; j++)
        for(k=kst; k<=ked; k++)
          val[i][j][k]=val[i+iof][j+jof][k+kof];
      }
    }
  }
}
