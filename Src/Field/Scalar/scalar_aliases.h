#ifndef ALIASES_H
#define ALIASES_H

#include "../../Boundary/bndcnd.h"
#include "../../Domain/domain.h"

/////////////
//         //
//  Shape  //
//         //
/////////////
struct Shape {
  Shape(const int ni, const int nj, const int nk, 
        const int oi, const int oj, const int ok,
        const int si, const int sj, const int sk, 
        const int ei, const int ej, const int ek,
        BndCnd * b, const Domain * d) :
   n_i(ni), n_j(nj), n_k(nk), o_i(oi), o_j(oj), o_k(ok),
   s_i(si), s_j(sj), s_k(sk), e_i(ei), e_j(ej), e_k(ek),
   dm(d)
   {bc=b;}
  const int n_i,n_j,n_k,  o_i,o_j,o_k,  s_i,s_j,s_k,  e_i,e_j,e_k;
  BndCnd       * bc;
  const Domain * dm;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: scalar_aliases.h,v 1.5 2008/11/17 19:23:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
