#include "vector.h"

/******************************************************************************/
void Vector::bnd_insert ( const Dir d, real *** cp ) {

  // i-min
  if( (d == Dir::imin()) && (dom->coord(Comp::i()) == 0) )
    for_m(m) {
      if(m==Comp::u()){
        for_mjk(m, j, k){
          vec[m][si(m)-1][j][k] = cp[~m][j][k];
        }
      } else {
        for_mjk(m, j, k){
          vec[m][si(m)  ][j][k] = cp[~m][j][k];
        }
      }
    }

  // i-max
  if( (d == Dir::imax()) && (dom->coord(Comp::i()) == dom->dim(Comp::i())-1) )
    for_m(m) {
      if(m==Comp::u()){
        for_mjk(m, j, k) 
          vec[m][ei(m)+1][j][k] = cp[~m][j][k];
      } else {
        for_mjk(m, j, k) 
          vec[m][ei(m)  ][j][k] = cp[~m][j][k];
      }
    }

  // j-min
  if( (d == Dir::jmin()) && (dom->coord(Comp::j()) == 0) )
    for_m(m) {
      if(m==Comp::v()){
        for_mik(m, i, k)
          vec[m][i][sj(m)-1][k] = cp[~m][i][k];
      } else {
        for_mik(m, i, k)
          vec[m][i][sj(m)  ][k] = cp[~m][i][k];
      }
    }

  // j-max
  if( (d == Dir::jmax()) && (dom->coord(Comp::j()) == dom->dim(Comp::j())-1) )
    for_m(m) {
      if(m==Comp::v()){
        for_mik(m, i, k)
          vec[m][i][ej(m)+1][k] = cp[~m][i][k];
      } else {
        for_mik(m, i, k)
          vec[m][i][ej(m)  ][k] = cp[~m][i][k];
      }
    }

  // k-min
  if( (d == Dir::kmin()) && (dom->coord(Comp::k()) == 0) )
    for_m(m) {
      if(m==Comp::w()){
        for_mij(m, i, j)
          vec[m][i][j][sk(m)-1] = cp[~m][i][j];
      } else {
        for_mij(m, i, j)
          vec[m][i][j][sk(m)  ] = cp[~m][i][j];
      }
    }

  // k-max
  if( (d == Dir::kmax()) && (dom->coord(Comp::k()) == dom->dim(Comp::k())-1) )
    for_m(m) {
      if(m==Comp::w()){
        for_mij(m, i, j)
          vec[m][i][j][ek(m)+1] = cp[~m][i][j];
      } else {
        for_mij(m, i, j)
          vec[m][i][j][ek(m)  ] = cp[~m][i][j];
      }
    }

  /* free the memory used during this process */
  dealloc3d( &cp );

  exchange();
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_bnd_insert.cpp,v 1.2 2016/03/15 15:38:24 sato Exp $'/
+-----------------------------------------------------------------------------*/
