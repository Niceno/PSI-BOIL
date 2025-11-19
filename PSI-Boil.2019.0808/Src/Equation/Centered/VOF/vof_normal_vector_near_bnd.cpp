#include "vof.h"

/******************************************************************************/
void VOF::normal_vector_near_bnd(const Scalar & val,const NormMethod & nm) {
/***************************************************************************//**
*  \brief normal vector for cells adjacent to a wall or an immersed boundary
*         obtained using extrapolated volume fractions
*******************************************************************************/
  
  /* in this case, no extraction of alpha is performed */
  Scalar * nalpha_ptr = &stmp;

  void(VOF::*norm_kernel) (real &, real &, real &, real &,
                           const int, const int, const int,
                           const Scalar &);
#if 1
  if       (nm==NormMethod::Mixed()) {
    norm_kernel = &VOF::norm_mixed_kernel;
  } else if(nm==NormMethod::Young()) {
    norm_kernel = &VOF::norm_young_kernel;
  } else if(nm==NormMethod::CC()) {
    norm_kernel = &VOF::norm_cc_kernel;
  } else if(nm==NormMethod::ElviraXZ()) {
    norm_kernel = &VOF::norm_elvira_kernel;
  } else if(nm==NormMethod::ElviraXY()) {
    norm_kernel = &VOF::norm_elvira_kernel;
  } else if(nm==NormMethod::ElviraYZ()) {
    norm_kernel = &VOF::norm_elvira_kernel;
  } else {
    /* default */
    boil::aout<<"VOF::normvector_near_wall: Normal vector calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }
#endif


  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::wall() ) {

      /*-------+
      |  Wall  |
      +-------*/

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      if(d != Dir::undefined()) {

        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1;
            //norm_kernel(nx[ii][j][k], ny[ii][j][k], nz[ii][j][k], ii,j ,k , val);
            (this->*norm_kernel)(nx[ii][j][k], ny[ii][j][k], nz[ii][j][k], 
                                 (*nalpha_ptr)[ii][j][k],
                                 ii,j ,k , val);
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;
            //norm_kernel(nx[ii][j][k], ny[ii][j][k], nz[ii][j][k], ii,j ,k , val);
            (this->*norm_kernel)(nx[ii][j][k], ny[ii][j][k], nz[ii][j][k],
                                 (*nalpha_ptr)[ii][j][k],
                                 ii,j ,k , val);
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
            //norm_kernel(nx[i][jj][k], ny[i][jj][k], nz[i][jj][k], i ,jj,k , val);
            (this->*norm_kernel)(nx[i][jj][k], ny[i][jj][k], nz[i][jj][k],
                                 (*nalpha_ptr)[i][jj][k],
                                 i ,jj,k , val);
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
            //norm_kernel(nx[i][jj][k], ny[i][jj][k], nz[i][jj][k], i ,jj,k , val);
            (this->*norm_kernel)(nx[i][jj][k], ny[i][jj][k], nz[i][jj][k],
                                 (*nalpha_ptr)[i][jj][k],
                                 
                                 i ,jj,k , val);
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
            //norm_kernel(nx[i][j][kk], ny[i][j][kk], nz[i][j][kk], i ,j ,kk, val);
            (this->*norm_kernel)(nx[i][j][kk], ny[i][j][kk], nz[i][j][kk],
                                 (*nalpha_ptr)[i][j][kk],
                                 
                                 i ,j ,kk, val);
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
            //norm_kernel(nx[i][j][kk], ny[i][j][kk], nz[i][j][kk], i ,j ,kk, val);
            (this->*norm_kernel)(nx[i][j][kk], ny[i][j][kk], nz[i][j][kk],
                                 (*nalpha_ptr)[i][j][kk],
                                 i ,j ,kk, val);
          }
        }
      }
    }

  } /* bcs */

  /***************+
  | immersed body |
  +***************/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    // cell[i][j][k] is wall adjacent cells in fluid domain
    dom->ibody().ijk(cc,&i,&j,&k);

    // west is in solid domain
    if (dom->ibody().off(i-1,j,k)) {
      //norm_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
      (this->*norm_kernel)(nx[i][j][k], ny[i][j][k], nz[i][j][k],
                           (*nalpha_ptr)[i][j][k],
                           i,j,k, val);
    }

    // east
    if (dom->ibody().off(i+1,j,k)) {
      //norm_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
      (this->*norm_kernel)(nx[i][j][k], ny[i][j][k], nz[i][j][k],
                           (*nalpha_ptr)[i][j][k],
                           i,j,k, val);
    }

    // south
    if (dom->ibody().off(i,j-1,k)) {
      //norm_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
      (this->*norm_kernel)(nx[i][j][k], ny[i][j][k], nz[i][j][k],
                           (*nalpha_ptr)[i][j][k],
                           i,j,k, val);
    }

    // north
    if (dom->ibody().off(i,j+1,k)) {
      //norm_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
      (this->*norm_kernel)(nx[i][j][k], ny[i][j][k], nz[i][j][k],
                           (*nalpha_ptr)[i][j][k],
                           i,j,k, val);
    }

    // bottom
    if (dom->ibody().off(i,j,k-1)) {
      //norm_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
      (this->*norm_kernel)(nx[i][j][k], ny[i][j][k], nz[i][j][k],
                           (*nalpha_ptr)[i][j][k],
                           i,j,k, val);
    }

    // top
    if (dom->ibody().off(i,j,k+1)) {
      //norm_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
      (this->*norm_kernel)(nx[i][j][k], ny[i][j][k], nz[i][j][k],
                           (*nalpha_ptr)[i][j][k],
                           i,j,k, val);
    }
  }

#if 1
  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();
#endif

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  return;
}
