#include "vof.h"

/******************************************************************************/
void VOF::extract_alpha(Scalar & scp) {
/***************************************************************************//**
 \brief Calculate value of alpha in cells
    if there is no interface in the cell, unreal=yotta (=1e+24) is stored.
    plane: vm1*x + vm2*y + vm3*z = alpha
    output: nalpha
    EDIT: I think the alpha is valid in standardized space
*******************************************************************************/

  /* avoid singular cases */
  for_avijk(scp,i,j,k) {
    if(scp[i][j][k]==0.5) scp[i][j][k] += boil::pico;
  }

  /* calculate alpha value in the domain */
  /* assumes positive normal vector in normalized space */
  for_vijk(nalpha,i,j,k) {
    nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  }

  nalpha.exchange_all();

  return;
}

/******************************************************************************/
void VOF::extract_alpha_near_bnd(const Scalar & scp) {
/***************************************************************************//**
*  \brief alpha for cells adjacent to a wall or an immersed boundary
*******************************************************************************/

  for( int b=0; b<scp.bc().count(); b++ ) {

    if(scp.bc().type_decomp(b)) continue;

    if( scp.bc().type(b) == BndType::wall() ) {

      /*-------+
      |  Wall  |
      +-------*/

      //int iof=0, jof=0, kof=0;

      Dir d      = scp.bc().direction(b);

      if(d != Dir::undefined()) {

        if(d == Dir::imin()){
          for_vijk( scp.bc().at(b), i,j,k ){
            int ii=i+1;
            nalpha[ii][j][k] = alpha_val(scp[ii][j][k],
                                nx[ii][j][k],ny[ii][j][k],nz[ii][j][k]);
          }
        }
        if(d == Dir::imax()){
          for_vijk( scp.bc().at(b), i,j,k ){
            int ii=i-1;
            nalpha[ii][j][k] = alpha_val(scp[ii][j][k],
                                nx[ii][j][k],ny[ii][j][k],nz[ii][j][k]);
          }
        }
        if(d == Dir::jmin()){
          for_vijk( scp.bc().at(b), i,j,k ){
            int jj=j+1;
            nalpha[i][jj][k] = alpha_val(scp[i][jj][k],
                                nx[i][jj][k],ny[i][jj][k],nz[i][jj][k]);
          }
        }
        if(d == Dir::jmax()){
          for_vijk( scp.bc().at(b), i,j,k ){
            int jj=j-1;
            nalpha[i][jj][k] = alpha_val(scp[i][jj][k],
                                nx[i][jj][k],ny[i][jj][k],nz[i][jj][k]);
          }
        }
        if(d == Dir::kmin()){
          for_vijk( scp.bc().at(b), i,j,k ){
            int kk=k+1;
            nalpha[i][j][kk] = alpha_val(scp[i][j][kk],
                                nx[i][j][kk],ny[i][j][kk],nz[i][j][kk]);
          }
        }
        if(d == Dir::kmax()){
          for_vijk( scp.bc().at(b), i,j,k ){
            int kk=k-1;
            nalpha[i][j][kk] = alpha_val(scp[i][j][kk],
                                nx[i][j][kk],ny[i][j][kk],nz[i][j][kk]);
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
      nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                  nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }

    // east
    if (dom->ibody().off(i+1,j,k)) {
      nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                  nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }

    // south
    if (dom->ibody().off(i,j-1,k)) {
      nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                  nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }

    // north
    if (dom->ibody().off(i,j+1,k)) {
      nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                  nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }

    // bottom
    if (dom->ibody().off(i,j,k-1)) {
      nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                  nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }

    // top
    if (dom->ibody().off(i,j,k+1)) {
      nalpha[i][j][k] = alpha_val(scp[i][j][k],
                                  nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }
  }

  nalpha.exchange_all();

  return;
}

/***********************
 * ancillary function
 ***********************/
real VOF::alpha_val(const real c, const real nnx, const real nny, const real nnz) {

  /* degenerate case I */
  if(c<boil::pico||c-1.0>-boil::pico) {
    return boil::unreal;
  }

  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -nnx;
  real vn2 = -nny;
  real vn3 = -nnz;

  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real denom = vm1+vm2+vm3;
  /* degenerate case II */
  if(denom<boil::pico)
    return boil::unreal;

  real qa = 1.0/denom;
  vm1 *= qa;
  vm2 *= qa;
  vm3 *= qa;

  return denom*calc_alpha(c, vm1, vm2, vm3);
}

