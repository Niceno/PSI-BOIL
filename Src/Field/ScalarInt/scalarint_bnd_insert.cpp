#include "scalarint.h"

/******************************************************************************/
void ScalarInt::bnd_insert ( const Dir d, int **  cp ) {

  /*-----------------------------------+
  |  direction is either imin or imax  |
  +-----------------------------------*/
  if( (d == Dir::imin()) && (dom->coord(Comp::i()) == 0) ) 
    for_ajk(j, k) 
      val[si()-1][j][k] = cp[j][k];

  if( (d == Dir::imax()) && (dom->coord(Comp::i()) == dom->dim(Comp::i())-1) )
    for_ajk(j, k) 
      val[ei()+1][j][k] = cp[j][k];

  /*-----------------------------------+
  |  direction is either jmin or jmax  |
  +-----------------------------------*/
  if( (d == Dir::jmin()) && (dom->coord(Comp::j()) == 0) ) 
    for_aik(i, k) 
      val[i][sj()-1][k] = cp[i][k];

  if( (d == Dir::jmax()) && (dom->coord(Comp::j()) == dom->dim(Comp::j())-1) )
    for_aik(i, k) 
      val[i][ej()+1][k] = cp[i][k];

  /*-----------------------------------+
  |  direction is either kmin or kmax  |
  +-----------------------------------*/
  if( (d == Dir::kmin()) && (dom->coord(Comp::k()) == 0) ) 
    for_aik(i, j) 
      val[i][j][sk()-1] = cp[i][j];

  if( (d == Dir::kmax()) && (dom->coord(Comp::k()) == dom->dim(Comp::k())-1) )
    for_aik(i, j) 
      val[i][j][ek()+1] = cp[i][j];

  /* free the memory used during this process */
  dealloc2d( &cp );
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarint_bnd_insert.cpp,v 1.1 2015/05/05 14:36:01 sato Exp $'/
+-----------------------------------------------------------------------------*/
