#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::bnd_insert ( const Dir d, real **  cp ) {

  /*-----------------------------------+
  |  direction is either imin or imax  |
  +-----------------------------------*/
  if( (d == Dir::imin()) && (dom->coord(Comp::i()) == 0) )
    for_ajk(j, k){
      clr[si()-1][j][k] = cp[j][k];
      phi[si()-1][j][k] = cp[j][k];
    }

  if( (d == Dir::imax()) && (dom->coord(Comp::i()) == dom->dim(Comp::i())-1) )
    for_ajk(j, k){
      clr[ei()+1][j][k] = cp[j][k];
      phi[ei()+1][j][k] = cp[j][k];
    }

  /*-----------------------------------+
  |  direction is either jmin or jmax  |
  +-----------------------------------*/
  if( (d == Dir::jmin()) && (dom->coord(Comp::j()) == 0) ) 
    for_aik(i, k){
      clr[i][sj()-1][k] = cp[i][k];
      phi[i][sj()-1][k] = cp[i][k];
    }

  if( (d == Dir::jmax()) && (dom->coord(Comp::j()) == dom->dim(Comp::j())-1) )
    for_aik(i, k){
      clr[i][ej()+1][k] = cp[i][k];
      phi[i][ej()+1][k] = cp[i][k];
    }

  /*-----------------------------------+
  |  direction is either kmin or kmax  |
  +-----------------------------------*/
  if( (d == Dir::kmin()) && (dom->coord(Comp::k()) == 0) ) 
    for_aik(i, j){
      clr[i][j][sk()-1] = cp[i][j];
      phi[i][j][sk()-1] = cp[i][j];
    }

  if( (d == Dir::kmax()) && (dom->coord(Comp::k()) == dom->dim(Comp::k())-1) )
    for_aik(i, j){
      clr[i][j][ek()+1] = cp[i][j];
      phi[i][j][ek()+1] = cp[i][j];
    }

  /* free the memory used during this process */
  dealloc2d( &cp );

  scheme.bdcond_f(clr);
  scheme.bdcond_i(clr);
  scheme.bdcond_j(clr);
  scheme.bdcond_k(clr);
  bdphiface(sxyz,Comp::i(),clr);
  bdphiface(sxyz,Comp::j(),clr);
  bdphiface(sxyz,Comp::k(),clr);
#ifdef IB
  ib_ext_scalar(clr);
  ib_bdcond(clr);
#endif
  //scheme.f.exchange_all();
  //scheme.sigx.exchange_all();
  //scheme.sigy.exchange_all();
  //scheme.sigz.exchange_all();
  //sxyz.exchange();

}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_bnd_insert.cpp,v 1.1 2015/01/05 17:10:52 sato Exp $'/
+-----------------------------------------------------------------------------*/
