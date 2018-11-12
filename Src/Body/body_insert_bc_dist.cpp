#include "body.h"
#include "../Field/Scalar/scalar.h"

/******************************************************************************/
void Body::insert_bc_dist(const Scalar * val) {
/***************************************************************************//**
*  \brief set boundary value for a distance function, etc.
*         dirichlet boundary condition is neglected.
*         scalar_exchange(_all) should take account of periodic condition.
*           1st: insert_bc_dist(phi);
*           2nd: phi.exchange_all();
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val->bc().count(); b++ ) {

    //std::cout<<val->bc().type(b)<<"\n";

    if( val->bc().type_decomp(b) ) continue;

    if( val->bc().type(b) == BndType::neumann()
      ||val->bc().type(b) == BndType::symmetry() ){

      int iof=0, jof=0, kof=0;

      Dir d      = val->bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d == Dir::ibody()){

      } else if(d != Dir::undefined()) {

        int ist, ied, jst, jed, kst, ked;
        ist=0;           jst=0;           kst=0;
	ied=val->ei()+1; jed=val->ej()+1; ked=val->ek()+1;

        if(d == Dir::imin()){ iof++; ist=val->si()-1; ied=val->si()-1;}
       	if(d == Dir::imax()){ iof--; ist=val->ei()+1; ied=val->ei()+1;}
        if(d == Dir::jmin()){ jof++; jst=val->sj()-1; jed=val->sj()-1;}
       	if(d == Dir::jmax()){ jof--; jst=val->ej()+1; jed=val->ej()+1;}
        if(d == Dir::kmin()){ kof++; kst=val->sk()-1; ked=val->sk()-1;}
       	if(d == Dir::kmax()){ kof--; kst=val->ek()+1; ked=val->ek()+1;}

        //for_vijk( val->bc().at(b), i,j,k )
	for(int i=ist; i<=ied; i++)
	  for(int j=jst; j<=jed; j++)
	    for(int k=kst; k<=ked; k++)
          (*val)[i][j][k]=(*val)[i+iof][j+jof][k+kof];
      }
    }

    if( val->bc().type(b) == BndType::wall()
      ||val->bc().type(b) == BndType::dirichlet()
      ||val->bc().type(b) == BndType::inlet()
      ||val->bc().type(b) == BndType::outlet()
      ||val->bc().type(b) == BndType::insert() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = val->bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d == Dir::ibody()){

      } else if(d != Dir::undefined()) {

        int ist, ied, jst, jed, kst, ked;
        ist=0;           jst=0;           kst=0;
	ied=val->ei()+1; jed=val->ej()+1; ked=val->ek()+1;

        if(d == Dir::imin()){ iof++; ist=val->si()-1; ied=val->si()-1;}
       	if(d == Dir::imax()){ iof--; ist=val->ei()+1; ied=val->ei()+1;}
        if(d == Dir::jmin()){ jof++; jst=val->sj()-1; jed=val->sj()-1;}
       	if(d == Dir::jmax()){ jof--; jst=val->ej()+1; jed=val->ej()+1;}
        if(d == Dir::kmin()){ kof++; kst=val->sk()-1; ked=val->sk()-1;}
       	if(d == Dir::kmax()){ kof--; kst=val->ek()+1; ked=val->ek()+1;}

        //for_vijk( val->bc().at(b), i,j,k )
	for(int i=ist; i<=ied; i++)
	  for(int j=jst; j<=jed; j++)
	    for(int k=kst; k<=ked; k++)
              (*val)[i][j][k]=(*val)[i+iof][j+jof][k+kof];
              //(*val)[i][j][k]=1.5*(*val)[i+iof][j+jof][k+kof]
              //               -0.5*(*val)[i+2*iof][j+2*jof][k+2*kof];
      }
    }

  }
}
