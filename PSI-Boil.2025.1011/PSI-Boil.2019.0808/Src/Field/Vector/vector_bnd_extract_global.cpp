#include "vector.h"

/******************************************************************************/
void Vector::bnd_extract_global( const Dir d, real **** cplane,
                                 const Vector * uvw_l ) const { 
  /* This function is used for velocity field without domain decomposition.
     uvw_l is velocity field with domain decomposition.
     Example:
       main: uvw_g.bnd_extract_global(Dir::imin(), &copy_planeVec, &uvw_l);
     where uvw_g and uvw_l are without and with decomposition, respectively. */ 

  /*-----------------------------------+
  |  direction is either imin or imax  |
  +-----------------------------------*/
  if( (d == Dir::imin()) || (d == Dir::imax()) ) {

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, uvw_l->nj()+1, uvw_l->nk()+1);

    if( d == Dir::imin() ) {
      for_m(m){
        //std::cout<<"bnd_extract:"<<m<<" "<<si(m)<<"\n";
        for_avmjk((*uvw_l), m, j, k) {
          int jj = uvw_l->dom->global_J(j);
          int kk = uvw_l->dom->global_K(k);
          (*cplane)[~m][j][k] = vec[m][si(m)][jj][kk];
        }
      }
    } else {
    for_m(m){
        //std::cout<<"bnd_extract:"<<m<<" "<<ei(m)<<"\n";
        for_avmjk((*uvw_l), m, j, k) {
          int jj = uvw_l->dom->global_J(j);
          int kk = uvw_l->dom->global_K(k);
          (*cplane)[~m][j][k] = vec[m][ei(m)][jj][kk];
        }
      }
    }
  }

  /*-----------------------------------+
  |  direction is either jmin or jmax  |
  +-----------------------------------*/
  if( (d == Dir::jmin()) || (d == Dir::jmax()) ) {

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, uvw_l->ni()+1, uvw_l->nk()+1);

    if( d == Dir::jmin() ) {
      for_m(m) {
        for_avmik((*uvw_l), m, i, k) {
          int ii = uvw_l->dom->global_I(i);
          int kk = uvw_l->dom->global_K(k);
          (*cplane)[~m][i][k] = vec[m][ii][sj(m)][kk];
        }
      }
    } else {
      for_m(m) {
        for_avmik((*uvw_l), m, i, k) {
          int ii = uvw_l->dom->global_I(i);
          int kk = uvw_l->dom->global_K(k);
          (*cplane)[~m][i][k] = vec[m][ii][ej(m)][kk];
        }       
      } 
    }
  }

  /*-----------------------------------+
  |  direction is either kmin or kmax  |
  +-----------------------------------*/
  if( (d == Dir::kmin()) || (d == Dir::kmax()) ) {

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, uvw_l->ni()+1, uvw_l->nj()+1);

    if( d == Dir::kmin() ) {
      for_m(m) {
        for_avmij((*uvw_l), m, i, j) {
          int ii = uvw_l->dom->global_I(i);
          int jj = uvw_l->dom->global_J(j);
          (*cplane)[~m][i][j] = vec[m][ii][jj][sk(m)];
        }
      }
    } else {
      for_m(m) {
        for_avmij((*uvw_l), m, i, j) {
          int ii = uvw_l->dom->global_I(i);
          int jj = uvw_l->dom->global_J(j);
          (*cplane)[~m][i][j] = vec[m][ii][jj][ek(m)];
        }
      }
    }
  }
}
