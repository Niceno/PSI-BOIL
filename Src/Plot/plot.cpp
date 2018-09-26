#include "plot.h"
#include "../Parallel/communicator.h"

/*-----------------+
|  global plotter  |
+-----------------*/
namespace boil {
  Plot * plot;
}

/*============================================================================*/
void Plot :: box_set(const Domain & d, const Box<int> & b) {
   
  /* set logical "boxed" to true */
  boxed=true;  

  /* copy the box */
  box = b;

  int g_vol_0 = 0;
  g_vol_0 = ( box.last (Comp::i()) - box.first(Comp::i()) + 1 )
          * ( box.last (Comp::j()) - box.first(Comp::j()) + 1 )
          * ( box.last (Comp::k()) - box.first(Comp::k()) + 1 );

  box.first( Comp::i(), d.local_i( box.first(Comp::i()) ) );
  box.last ( Comp::i(), d.local_i( box.last (Comp::i()) ) );
  box.first( Comp::j(), d.local_j( box.first(Comp::j()) ) );
  box.last ( Comp::j(), d.local_j( box.last (Comp::j()) ) );
  box.first( Comp::k(), d.local_k( box.first(Comp::k()) ) );
  box.last ( Comp::k(), d.local_k( box.last (Comp::k()) ) );

  for_m(m) {
    if( box.first(m) == -1 && box.last(m) == -1 ) {
      /* box is outside the domain */
      box.first(m, 0);
      box.last (m,-1);
      /* box encloses the subdoimain */
      if(m==Comp::i()) 
        if( b.first(m) <= d.global_I(1) && b.last(m) >= d.global_I(d.ni()-2) ) {
          box.first(m,1);
          box.last(m,d.ni()-2);
        }
      if(m==Comp::j()) 
        if( b.first(m) <= d.global_J(1) && b.last(m) >= d.global_J(d.nj()-2) ) {
          box.first(m,1);
          box.last(m,d.nj()-2);
        }
      if(m==Comp::k()) 
        if( b.first(m) <= d.global_K(1) && b.last(m) >= d.global_K(d.nk()-2) ) {
          box.first(m,1);
          box.last(m,d.nk()-2);
        }
    }
    else if( box.first(m) == -1 && box.last(m) != -1 ) {
      box.first(m, 1);
    }
    else if( box.first(m) != -1 && box.last(m) == -1 ) {
      if(m==Comp::i()) box.last(m, d.ni()-2);
      if(m==Comp::j()) box.last(m, d.nj()-2);
      if(m==Comp::k()) box.last(m, d.nk()-2);
    }
  }

  int g_vol_1 = 0;
  g_vol_1 = ( box.last (Comp::i()) - box.first(Comp::i()) + 1 )
          * ( box.last (Comp::j()) - box.first(Comp::j()) + 1 )
          * ( box.last (Comp::k()) - box.first(Comp::k()) + 1 );
  boil::cart.sum_int(&g_vol_1);

  assert( g_vol_1 == g_vol_0 );
}

/*============================================================================*/
void Plot :: box_release() {
  boxed=false;
}

/*-----------------------------------------------------------------------------+
 '$Id: plot.cpp,v 1.5 2011/11/03 19:24:28 niceno Exp $'/
+-----------------------------------------------------------------------------*/
