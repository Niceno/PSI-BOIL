#include "vector.h"
//#define DEBUG

/******************************************************************************/
void Vector::bnd_extract_pos( const Dir d, real **** cplane,
                              const Vector * uvw_l ) const { 
  /* This function is used to extract velocity field based on position defined.
     uvw_l is velocity field with domain decomposition.
     Example:
     main: uvw_o.bnd_extract_global(Dir::imin(), &copy_planeVec, &uvw_l);
     where uvw_o is the origin, and uvw_l is the destination.  */ 

  real *** oplane; // origine plane
  int  *** oplane_int; // origine plane
  int gi = dom->gi(); // global ni of origin
  int gj = dom->gj(); // global nj of origin
  int gk = dom->gk(); // global nk of origin

#ifdef DEBUG
  std::cout<<"bnd_extract_pos: gi= "<<gi<<" gj= "<<gj<<" gk= "<<gk<<"\n";
#endif

  /*-----------------------------------+
  |  direction is either imin or imax  |
  +-----------------------------------*/
  if( (d == Dir::imin()) || (d == Dir::imax()) ) {

    /* allocate memory for the origine plane */
    alloc3d( &oplane, 3, gj, gk);
    alloc3d( &oplane_int, 3, gj, gk);

    /* initialize */
    for(int m=0; m<3; m++)
      for(int j=0; j<gj; j++)
        for(int k=0; k<gk; k++) {
          oplane    [m][j][k]=0.0;
          oplane_int[m][j][k]=0;
        }

    if( (d == Dir::imin()) &&
        (dom->coord(Comp::i()) == 0) ) {
      //std::cout<<"imin:: rank= "<<boil::cart.iam()<<"\n";
      for_m(m){
        for_mjk( m, j, k) {
          int J = dom->global_J(j);
          int K = dom->global_K(k);
          oplane    [~m][J][K] = vec[m][si(m)][j][k];
          oplane_int[~m][J][K] = 1.0;
        }
      }
    }
    if( (d == Dir::imax()) && 
        (dom->coord(Comp::i()) == dom->dim(Comp::i())-1) ) {
      //std::cout<<"imax:: rank= "<<boil::cart.iam()<<"\n";
      for_m(m){
        //std::cout<<boil::cart.iam()<<" "<<m<<" "<<ei(m)<<"\n";
        for_mjk( m, j, k) {
          int J = dom->global_J(j);
          int K = dom->global_K(k);
          oplane    [~m][J][K] = vec[m][ei(m)][j][k];
          oplane_int[~m][J][K] = 1.0;
        }
      }
    }
    boil::cart.sum_real_n(**oplane, 3*gj*gk);
    boil::cart.sum_int_n (**oplane_int, 3*gj*gk);

    for(int m=0; m<3; m++)
      for(int j=0; j<gj; j++)
        for(int k=0; k<gk; k++) {
          if (oplane_int[m][j][k]!=0) {
            oplane[m][j][k] /= real(oplane_int[m][j][k]);
          }
        }
#if 0
    if (boil::cart.iam()==0) {
      for(int k=0; k<gk; k++) {
        std::cout<<k<<" "<<oplane[0][1][k]<<"\n";
      }
    }
#endif

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, uvw_l->nj()+1, uvw_l->nk()+1);

    for_m(m){
      real ys = yc_global(m,1);
      real ye = yc_global(m,dom->gj()-2);
      real zs = zc_global(m,1);
      real ze = zc_global(m,dom->gk()-2);
      if (m==Comp::j()) ye = yc_global(m,dom->gj()-1);
      if (m==Comp::k()) ze = zc_global(m,dom->gk()-1);
      for_vmjk((*uvw_l), m, j, k) {
        bool inside = true;
        real yy = uvw_l->yc(m,j);
        real zz = uvw_l->zc(m,k);
        if ( yy < ys-boil::pico) inside = false;
        if ( yy > ye+boil::pico) inside = false;
        if ( zz < zs-boil::pico) inside = false;
        if ( zz > ze+boil::pico) inside = false;
        if (inside) {
          int JJ = J(m,yy); 
          int KK = K(m,zz); 
          (*cplane)[~m][j][k] = oplane[~m][JJ][KK];
        } else {
          (*cplane)[~m][j][k] = 0.0;
        }
      }
    }
    dealloc3d( &oplane );
    dealloc3d( &oplane_int );
  }

  /*-----------------------------------+
  |  direction is either jmin or jmax  |
  +-----------------------------------*/
  if( (d == Dir::jmin()) || (d == Dir::jmax()) ) {

    /* allocate memory for the origine plane */
    alloc3d( &oplane, 3, gi, gk);
    alloc3d( &oplane_int, 3, gi, gk);

    /* initialize */
    for(int m=0; m<3; m++)
      for(int i=0; i<gi; i++)
        for(int k=0; k<gk; k++) {
          oplane    [m][i][k]=0.0;
          oplane_int[m][i][k]=0;
        }

    if( (d == Dir::jmin()) &&
        (dom->coord(Comp::j()) == 0) ) {
      for_m(m){
        for_mik( m, i, k) {
          int I = dom->global_I(i);
          int K = dom->global_K(k);
          oplane    [~m][I][K] = vec[m][i][sj(m)][k];
          oplane_int[~m][I][K] = 1.0;
        }
      }
    }
    if( (d == Dir::jmax()) && 
        (dom->coord(Comp::j()) == dom->dim(Comp::j())-1) ) {
      for_m(m){
        for_mik( m, i, k) {
          int I = dom->global_I(i);
          int K = dom->global_K(k);
          oplane    [~m][I][K] = vec[m][i][ej(m)][k];
          oplane_int[~m][I][K] = 1.0;
        }
      }
    }
    boil::cart.sum_real_n(**oplane, 3*gi*gk);
    boil::cart.sum_int_n (**oplane_int, 3*gi*gk);

    for(int m=0; m<3; m++)
      for(int i=0; i<gi; i++)
        for(int k=0; k<gk; k++) {
          if (oplane_int[m][i][k]!=0) {
            oplane[m][i][k] /= real(oplane_int[m][i][k]);
          }
        }

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, uvw_l->ni()+1, uvw_l->nk()+1);

    for_m(m){
      real xs = xc_global(m,1);
      real xe = xc_global(m,dom->gi()-2);
      real zs = zc_global(m,1);
      real ze = zc_global(m,dom->gk()-2);
      if (m==Comp::i()) xe = xc_global(m,dom->gi()-1);
      if (m==Comp::k()) ze = zc_global(m,dom->gk()-1);
      for_vmik((*uvw_l), m, i, k) {
        bool inside = true;
        real xx = uvw_l->xc(m,i);
        real zz = uvw_l->zc(m,k);
        if ( xx < xs-boil::pico) inside = false;
        if ( xx > xe+boil::pico) inside = false;
        if ( zz < zs-boil::pico) inside = false;
        if ( zz > ze+boil::pico) inside = false;
        if (inside) {
          int II = I(m,xx); 
          int KK = K(m,zz); 
          (*cplane)[~m][i][k] = oplane[~m][II][KK];
        } else {
          (*cplane)[~m][i][k] = 0.0;
        }
      }
    }
    dealloc3d( &oplane );
    dealloc3d( &oplane_int );
  }

  /*-----------------------------------+
  |  direction is either kmin or kmax  |
  +-----------------------------------*/
  if( (d == Dir::kmin()) || (d == Dir::kmax()) ) {

    /* allocate memory for the origine plane */
    alloc3d( &oplane, 3, gi, gj);
    alloc3d( &oplane_int, 3, gi, gj);

    /* initialize */
    for(int m=0; m<3; m++)
      for(int i=0; i<gi; i++)
        for(int j=0; j<gj; j++) {
          oplane    [m][i][j]=0.0;
          oplane_int[m][i][j]=0;
        }

    if( (d == Dir::kmin()) &&
        (dom->coord(Comp::k()) == 0) ) {
      for_m(m){
        for_mij( m, i, j) {
          int I = dom->global_I(i);
          int J = dom->global_J(j);
          oplane    [~m][I][J] = vec[m][i][j][sk(m)];
          oplane_int[~m][I][J] = 1.0;
        }
      }
    }
    if( (d == Dir::kmax()) && 
        (dom->coord(Comp::k()) == dom->dim(Comp::k())-1) ) {
      for_m(m){
        for_mij( m, i, j) {
          int I = dom->global_I(i);
          int J = dom->global_J(j);
          oplane    [~m][I][J] = vec[m][i][j][ek(m)];
          oplane_int[~m][I][J] = 1.0;
        }
      }
    }
    boil::cart.sum_real_n(**oplane, 3*gi*gj);
    boil::cart.sum_int_n (**oplane_int, 3*gi*gj);

    for(int m=0; m<3; m++)
      for(int i=0; i<gi; i++)
        for(int j=0; j<gj; j++)
          if (oplane_int[m][i][j]!=0) {
            oplane[m][i][j] /= real(oplane_int[m][i][j]);
          }

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, uvw_l->ni()+1, uvw_l->nj()+1);

    for_m(m){
      real xs = xc_global(m,1);
      real xe = xc_global(m,dom->gi()-2);
      real ys = yc_global(m,1);
      real ye = yc_global(m,dom->gj()-2);
      if (m==Comp::i()) xe = xc_global(m,dom->gi()-1);
      if (m==Comp::j()) ye = yc_global(m,dom->gj()-1);
      for_vmij((*uvw_l), m, i, j) {
        bool inside = true;
        real xx = uvw_l->xc(m,i);
        real yy = uvw_l->yc(m,j);
        if ( xx < xs-boil::pico) inside = false;
        if ( xx > xe+boil::pico) inside = false;
        if ( yy < ys-boil::pico) inside = false;
        if ( yy > ye+boil::pico) inside = false;
        if (inside) {
          int II = I(m,xx); 
          int JJ = J(m,yy); 
          (*cplane)[~m][i][j] = oplane[~m][II][JJ];
        } else {
          (*cplane)[~m][i][j] = 0.0;
        }
      }
    }
    dealloc3d( &oplane );
    dealloc3d( &oplane_int );
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_bnd_extract_pos.cpp,v 1.1 2016/03/15 15:39:26 sato Exp $'/
+-----------------------------------------------------------------------------*/
