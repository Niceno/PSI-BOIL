#include "scalar.h"

/******************************************************************************/
void Scalar::exchange_all_avg(const int dir) const {
/*------------------------------------------------------------------+
|  this functions averages velocities at the processor boundaries.  |
|  it is important for parallel version and ONLY for velocities.    |
+------------------------------------------------------------------*/
	
  if( bc().count() == 0 ) {
    OMS(Warning: exchanging a variable without boundary conditions);
  }

  /*-------------------------------+
  |                                |
  |  buffers for parallel version  |
  |                                |
  +-------------------------------*/
  real * sbuff_s, * sbuff_e, * rbuff_s, * rbuff_e;

  par_request req_s1,req_s2,req_r1,req_r2;
  par_status  ps;

  /* allocate memory for buffers */
  const int n = boil::maxi(ni(), nj(), nk())+1;
  assert(n > 0);

  sbuff_s = new real [n*n];
  sbuff_e = new real [n*n];
  rbuff_s = new real [n*n];
  rbuff_e = new real [n*n];

  /*----------------+
  |  I - direction  |
  +----------------*/
  if((dir == -1 || dir == 0) && o_x == 1) {

    /* decomposed in I direction */
    if( dom->dim(Comp::i()) > 1 ) {
      for_ajk(j,k) {
        int l = k*nj()+j;
        sbuff_e[l] = val[e_x][j][k];   // buffer i end
        sbuff_s[l] = val[s_x][j][k];   // buffer i start
        rbuff_e[l] = val[e_x][j][k];
        rbuff_s[l] = val[s_x][j][k];
      }
      /* send last and receive first */
      boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], nj()*nk(),
                          par_real, dom->neighbour(Dir::imax()),
                                    dom->neighbour(Dir::imin()), Tag(0));

      /* send first and receive last */
      boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], nj()*nk(),
                          par_real, dom->neighbour(Dir::imin()),
                                    dom->neighbour(Dir::imax()), Tag(1));
      for_ajk(j,k) {
        int l = k*nj()+j;
        val[e_x][j][k] = 0.5*(sbuff_e[l]+rbuff_e[l]);   // buffer i end
        val[s_x][j][k] = 0.5*(sbuff_s[l]+rbuff_s[l]);   // buffer i start
      }
    } else if ( dom->dim(Comp::i()) == 1 ) {
      if( bc().type(Dir::imin(), BndType::periodic()) &&
          bc().type(Dir::imax(), BndType::periodic()) ) {
        for_ajk(j,k) {
          int l = k*nj()+j;
          sbuff_s[l] = val[s_x][j][k];   // buffer i start
          sbuff_e[l] = val[e_x][j][k];   // buffer i end
        }
        for_ajk(j,k) {
          int l = k*nj()+j;
          val[s_x][j][k] = 0.5*(sbuff_s[l]+sbuff_e[l]);
          val[e_x][j][k] = 0.5*(sbuff_s[l]+sbuff_e[l]);
        }
      }
    }
  }
  
  /*----------------+
  |  J - direction  |
  +----------------*/
  if((dir == -1 || dir == 1) && o_y == 1) {

    /* decomposed in J direction */
    if( dom->dim(Comp::j()) > 1 ) {
      for_aik(i,k) {
        int l = k*ni()+i;
        sbuff_e[l] = val[i][e_y][k];   // buffer j end
        sbuff_s[l] = val[i][s_y][k];   // buffer j start
        rbuff_e[l] = val[i][e_y][k];
        rbuff_s[l] = val[i][s_y][k];
      }
      /* send last and receive first */
      boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], ni()*nk(),
                          par_real, dom->neighbour(Dir::jmax()),
                                    dom->neighbour(Dir::jmin()), Tag(2));

      /* send first and receive last */
      boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], ni()*nk(),
                          par_real, dom->neighbour(Dir::jmin()),
                                    dom->neighbour(Dir::jmax()), Tag(3));
      for_aik(i,k) {
        int l = k*ni()+i;
        val[i][e_y][k] = 0.5*(sbuff_e[l]+rbuff_e[l]);
        val[i][s_y][k] = 0.5*(sbuff_s[l]+rbuff_s[l]);
      }
    } else if ( dom->dim(Comp::j()) == 1 ) {
      if( bc().type(Dir::jmin(), BndType::periodic()) &&
          bc().type(Dir::jmax(), BndType::periodic()) ) {
        for_aik(i,k) {
          int l = k*ni()+i;
          sbuff_s[l] = val[i][s_y][k];   // buffer j start
          sbuff_e[l] = val[i][e_y][k];   // buffer j end
        }
        for_aik(i,k) {
          int l = k*ni()+i;
          val[i][s_y][k] = 0.5*(sbuff_s[l]+sbuff_e[l]);
          val[i][e_y][k] = 0.5*(sbuff_s[l]+sbuff_e[l]);
        }
      }
    }
  }
  
  /*----------------+
  |  K - direction  |
  +----------------*/
  if((dir == -1 || dir == 2) && o_z == 1) {

    /* decomposed in K direction */
    if( dom->dim(Comp::k()) > 1 ) {
      for_aij(i,j) {
        int l = j*ni()+i;
        sbuff_e[l] = val[i][j][e_z];   // buffer k end
        sbuff_s[l] = val[i][j][s_z];   // buffer k start
        rbuff_e[l] = val[i][j][e_z];
        rbuff_s[l] = val[i][j][s_z];
      }
      /* send last and receive first */
      boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], ni()*nj(),
                          par_real, dom->neighbour(Dir::kmax()),
                                    dom->neighbour(Dir::kmin()), Tag(4));

      /* send first and receive last */
      boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], ni()*nj(),
                          par_real, dom->neighbour(Dir::kmin()),
                                    dom->neighbour(Dir::kmax()), Tag(5));
      for_aij(i,j) {
        int l = j*ni()+i;
        val[i][j][e_z] = 0.5*(sbuff_e[l]+rbuff_e[l]);   // buffer k end
        val[i][j][s_z] = 0.5*(sbuff_s[l]+rbuff_s[l]);   // buffer k start
      }
    } else if( dom->dim(Comp::k()) == 1 ) {
      if( bc().type(Dir::kmin(), BndType::periodic()) &&
          bc().type(Dir::kmax(), BndType::periodic()) ) {
        for_aij(i,j) {
          int l = j*ni()+i;
          sbuff_s[l] = val[i][j][s_z];   // buffer k start
          sbuff_e[l] = val[i][j][e_z];   // buffer k end
        }
        for_aij(i,j) {
          int l = j*ni()+i;
          val[i][j][s_z] = 0.5*(sbuff_s[l]+sbuff_e[l]);
          val[i][j][e_z] = 0.5*(sbuff_s[l]+sbuff_e[l]);
        }
      }
    }
  }
  
  delete [] sbuff_s;
  delete [] sbuff_e;
  delete [] rbuff_s;
  delete [] rbuff_e;
}

/*-----------------------------------------------------------------------------+
 '$Id: scalar_exchange_all_avg.cpp,v 1.5 2017/04/12 15:50:58 sato Exp $'/
+-----------------------------------------------------------------------------*/
