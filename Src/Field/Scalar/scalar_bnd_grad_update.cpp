#include "scalar.h"

/******************************************************************************/
void Scalar::bnd_grad_update(const Comp & com) {
/*-----------------------------------------------------------------------------+
|  boundary condition for dval/dx, dval/dy, dval/dz                            |
+-----------------------------------------------------------------------------*/

  boil::oout<<"Possible underdevelopment due to changes in buffer treatment ("
             <<"cutoff introduction). Review this function before using. "
             <<"Exiting."<<boil::endl;
  exit(0);

  if( bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0){
  } else if( bc().type_here(Dir::imin(),BndType::symmetry())
          //|| bc().type_here(Dir::imin(),BndType::wall())
          || bc().type_here(Dir::imin(),BndType::neumann())) {
    int i=si();
    if (com ==Comp::i()) {      // val=dphidx
      for_jk(j,k)
        //val[i-1][j][k]=0.0;
        val[i-1][j][k]=-val[i][j][k];
    } else {
      for_jk(j,k)
        val[i-1][j][k]=val[i][j][k];
    }
  } else {
    int i=si();
    if (com ==Comp::i()) {      // val=dphidx
      for_jk(j,k)
        val[i-1][j][k]=val[i][j][k];
    } else {
      for_jk(j,k)
        val[i-1][j][k]=val[i][j][k];
    }
  }

  if( bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) !=dom->dim(Comp::i())-1){
  } else if( bc().type_here(Dir::imax(),BndType::symmetry())
          //|| bc().type_here(Dir::imax(),BndType::wall())
          || bc().type_here(Dir::imax(),BndType::neumann())){
    int i=ei();
    if (com==Comp::i()) {       // val=dphidx
      for_jk(j,k)
        //val[i+1][j][k]=0.0;
        val[i+1][j][k]=-val[i][j][k];
    } else {
      for_jk(j,k)
        val[i+1][j][k]=val[i][j][k];
    }
  } else {
    int i=ei();
    if (com==Comp::i()) {       // val=dphidx
      for_jk(j,k)
        val[i+1][j][k]=val[i][j][k];
    } else {
      for_jk(j,k)
        val[i+1][j][k]=val[i][j][k];
    }
  }

  if( bc().type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0){
  } else if( bc().type_here(Dir::jmin(),BndType::symmetry())
          //|| bc().type_here(Dir::jmin(),BndType::wall())
          || bc().type_here(Dir::jmin(),BndType::neumann())) {
    int j=sj();
    if (com==Comp::j()) {       // val=dphidy
      for_ik(i,k)
        //val[i][j-1][k]=0.0;
        val[i][j-1][k]=-val[i][j][k];
    } else {
      for_ik(i,k)
        val[i][j-1][k]=val[i][j][k];
    }
  } else {
    int j=sj();
    if (com==Comp::j()) {       // val=dphidy
      for_ik(i,k)
        val[i][j-1][k]=val[i][j][k];
    } else {
      for_ik(i,k)
        val[i][j-1][k]=val[i][j][k];
    }
  }

  if( bc().type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) !=dom->dim(Comp::j())-1){
  } else if( bc().type_here(Dir::jmax(),BndType::symmetry())
          //|| bc().type_here(Dir::jmax(),BndType::wall())
          || bc().type_here(Dir::jmax(),BndType::neumann())) {
    int j=ej();
    if (com==Comp::j()) {       // val=dphidy
      for_ik(i,k)
        //val[i][j+1][k]=0.0;
        val[i][j+1][k]=-val[i][j][k];
    } else {
      for_ik(i,k)
        val[i][j+1][k]=val[i][j][k];
    }
  } else {
    int j=ej();
    if (com==Comp::j()) {       // val=dphidy
      for_ik(i,k)
        val[i][j+1][k]=val[i][j][k];
    } else {
      for_ik(i,k)
        val[i][j+1][k]=val[i][j][k];
    }
  }

  if( bc().type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0){
  } else if( bc().type_here(Dir::kmin(),BndType::symmetry())
          //|| bc().type_here(Dir::kmin(),BndType::wall())
          || bc().type_here(Dir::kmin(),BndType::neumann())) {
    int k=sk();
    if (com==Comp::k()) {       // val=dphidz
      for_ij(i,j)
        //val[i][j][k-1]=0.0;
        val[i][j][k-1]=-val[i][j][k];
    } else {
      for_ij(i,j)
        val[i][j][k-1]=val[i][j][k];
    }
  } else {
    int k=sk();
    if (com==Comp::k()) {       // val=dphidz
      for_ij(i,j)
        val[i][j][k-1]=val[i][j][k];
    } else {
      for_ij(i,j)
        val[i][j][k-1]=val[i][j][k];
    }
  }

  if( bc().type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1){
  } else if( bc().type_here(Dir::kmax(),BndType::symmetry())
          //|| bc().type_here(Dir::kmax(),BndType::wall())
          || bc().type_here(Dir::kmax(),BndType::neumann())){
    int k=ek();
    if (com==Comp::k()) {       // val=dphidz
      for_ij(i,j)
        //val[i][j][k+1]=0.0;
        val[i][j][k+1]=-val[i][j][k];
    } else {
      for_ij(i,j)
        val[i][j][k+1]=val[i][j][k];
    }
  } else {
    int k=ek();
    if (com==Comp::k()) {       // val=dphidz
      for_ij(i,j)
        val[i][j][k+1]=val[i][j][k];
    } else {
      for_ij(i,j)
        val[i][j][k+1]=val[i][j][k];
    }
  }

}
