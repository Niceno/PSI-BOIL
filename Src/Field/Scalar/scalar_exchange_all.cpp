#include "scalar.h"

/* uncomment the line below for blocking (old) send-receive */
//#define SENDRECV

/******************************************************************************/
void Scalar::exchange_all(const int dir) const {
	
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

  /*-------------------------------+
  |                                |
  |  browse through buffer layers  |
  |                                |
  +-------------------------------*/
  for(int b=boil::BW-1; b>=0; b--) {

  /*----------------+
  |  I - direction  |
  +----------------*/
  if(dir == -1 || dir == 0) {

    /* not decomposed in I direction */
    if( dom->dim(Comp::i()) == 1 ) {
      if( bc().type(Dir::imin(), BndType::periodic()) && 
          bc().type(Dir::imax(), BndType::periodic()) ) {
        int ib_min = bc().index(Dir::imin(), BndType::periodic());
        int ib_max = bc().index(Dir::imax(), BndType::periodic());
        real v_min = bc().value(ib_min);
        real v_max = bc().value(ib_max);
        for_ajk(j,k) {
          val[e_x + 1 + b][j][k] = val[s_x + o_x + b][j][k] + v_max;
          val[s_x - 1 - b][j][k] = val[e_x - o_x - b][j][k] + v_min;
        }
      }
    }
    /* decomposed */
    else {
#ifdef SENDRECV
      for_ajk(j,k) {
        int l = k*nj()+j;
        assert(l < n*n);
        sbuff_e[l] = val[e_x - o_x - b][j][k];   // buffer i end
        sbuff_s[l] = val[s_x + o_x + b][j][k];   // buffer i start
        rbuff_e[l] = val[e_x +  1  + b][j][k];
        rbuff_s[l] = val[s_x -  1  - b][j][k];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], nj()*nk(),
                            par_real, dom->neighbour(Dir::imax()),
                                      dom->neighbour(Dir::imin()), Tag(0));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], nj()*nk(),
                            par_real, dom->neighbour(Dir::imin()),
                                      dom->neighbour(Dir::imax()), Tag(1));
#else
      for_ajk(j,k) {
        int l = k*nj()+j;
        assert(l < n*n);
        rbuff_e[l] = val[e_x +  1 + b][j][k];
        rbuff_s[l] = val[s_x -  1 - b][j][k];
      }

      if( dom->neighbour(Dir::imin()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_s[0], nj()*nk(), par_real, 
                           dom->neighbour(Dir::imin()), Tag(0), & req_r1 );
      }
      if( dom->neighbour(Dir::imax()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_e[0], nj()*nk(), par_real, 
                           dom->neighbour(Dir::imax()), Tag(1), & req_r2 );
      }

      for_ajk(j,k) {
        int l = k*nj()+j;
        assert(l < n*n);
        sbuff_e[l] = val[e_x - o_x - b][j][k];   // buffer i end
        sbuff_s[l] = val[s_x + o_x + b][j][k];   // buffer i start
      }

      if( dom->neighbour(Dir::imax()) != par_proc_null ) {
        boil::cart.isend( &sbuff_e[0], nj()*nk(), par_real, 
                           dom->neighbour(Dir::imax()), Tag(0), & req_s1 );
      }
      if( dom->neighbour(Dir::imin()) != par_proc_null ) {
        boil::cart.isend( &sbuff_s[0], nj()*nk(), par_real, 
                           dom->neighbour(Dir::imin()), Tag(1), & req_s2 );
      }

      if( dom->neighbour(Dir::imax()) != par_proc_null )
        boil::cart.wait( & req_s1, & ps );
      if( dom->neighbour(Dir::imin()) != par_proc_null )
        boil::cart.wait( & req_s2, & ps );
      if( dom->neighbour(Dir::imin()) != par_proc_null )
        boil::cart.wait( & req_r1, & ps );
      if( dom->neighbour(Dir::imax()) != par_proc_null )
        boil::cart.wait( & req_r2, & ps );
#endif
      real v_min=0.0;
      real v_max=0.0;
      if(bc().type_here(Dir::imin(), BndType::periodic())){
        int ib_min = bc().index(Dir::imin(), BndType::periodic());
        v_min = bc().value(ib_min);
      }
      if(bc().type_here(Dir::imax(), BndType::periodic())){
        int ib_max = bc().index(Dir::imax(), BndType::periodic());
        v_max = bc().value(ib_max);
      }
      for_ajk(j,k) {
        int l = k*nj()+j;
        assert(l < n*n);
        val[e_x + 1 + b][j][k] = rbuff_e[l] + v_max;   // buffer i end
        val[s_x - 1 - b][j][k] = rbuff_s[l] + v_min;   // buffer i start
      }
    }
  }

  /*----------------+
  |  J - direction  |
  +----------------*/
  if(dir == -1 || dir == 1) {

    /* not decomposed in J direction */
    if( dom->dim(Comp::j()) == 1 ) {
      if( bc().type(Dir::jmin(), BndType::periodic()) && 
          bc().type(Dir::jmax(), BndType::periodic()) ) {
        int jb_min = bc().index(Dir::jmin(), BndType::periodic());
        int jb_max = bc().index(Dir::jmax(), BndType::periodic());
        real v_min = bc().value(jb_min);
        real v_max = bc().value(jb_max);
        for_aik(i,k) {
          val[i][e_y + 1 + b][k] = val[i][s_y + o_y + b][k] + v_max;
          val[i][s_y - 1 - b][k] = val[i][e_y - o_y - b][k] + v_min;
        }
      }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_aik(i,k) {
        int l = k*ni()+i;
        assert(l < n*n);
        sbuff_e[l] = val[i][e_y - o_y - b][k];   // buffer j end
        sbuff_s[l] = val[i][s_y + o_y + b][k];   // buffer j start
        rbuff_e[l] = val[i][e_y +  1  + b][k];
        rbuff_s[l] = val[i][s_y -  1  - b][k];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], ni()*nk(),
                            par_real, dom->neighbour(Dir::jmax()),
                                      dom->neighbour(Dir::jmin()), Tag(2));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], ni()*nk(),
                            par_real, dom->neighbour(Dir::jmin()),
                                      dom->neighbour(Dir::jmax()), Tag(3));
#else
      for_aik(i,k) {
        int l = k*ni()+i;
        assert(l < n*n);
        rbuff_e[l] = val[i][e_y +  1 + b][k];
        rbuff_s[l] = val[i][s_y -  1 - b][k];
      }

      if( dom->neighbour(Dir::jmin()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_s[0], ni()*nk(), par_real, 
                           dom->neighbour(Dir::jmin()), Tag(2), & req_r1 );
      }
      if( dom->neighbour(Dir::jmax()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_e[0], ni()*nk(), par_real, 
                           dom->neighbour(Dir::jmax()), Tag(3), & req_r2 );
      }

      for_aik(i,k) {
        int l = k*ni()+i;
        assert(l < n*n);
        sbuff_e[l] = val[i][e_y - o_y - b][k];   // buffer j end
        sbuff_s[l] = val[i][s_y + o_y + b][k];   // buffer j start
      }

      if( dom->neighbour(Dir::jmax()) != par_proc_null ) {
        boil::cart.isend( &sbuff_e[0], ni()*nk(), par_real, 
                           dom->neighbour(Dir::jmax()), Tag(2), & req_s1 );
      }
      if( dom->neighbour(Dir::jmin()) != par_proc_null ) {
        boil::cart.isend( &sbuff_s[0], ni()*nk(), par_real, 
                           dom->neighbour(Dir::jmin()), Tag(3), & req_s2 );
      }

      if( dom->neighbour(Dir::jmax()) != par_proc_null )
        boil::cart.wait( & req_s1, & ps );
      if( dom->neighbour(Dir::jmin()) != par_proc_null )
        boil::cart.wait( & req_s2, & ps );
      if( dom->neighbour(Dir::jmin()) != par_proc_null )
        boil::cart.wait( & req_r1, & ps );
      if( dom->neighbour(Dir::jmax()) != par_proc_null )
        boil::cart.wait( & req_r2, & ps );
#endif
      real v_min=0.0;
      real v_max=0.0;
      if(bc().type_here(Dir::jmin(), BndType::periodic())){
        int jb_min = bc().index(Dir::jmin(), BndType::periodic());
        v_min = bc().value(jb_min);
      }
      if(bc().type_here(Dir::jmax(), BndType::periodic())){
        int jb_max = bc().index(Dir::jmax(), BndType::periodic());
        v_max = bc().value(jb_max);
      }
      for_aik(i,k) {
        int l = k*ni()+i;
        assert(l < n*n);
        val[i][e_y + 1 + b][k] = rbuff_e[l] + v_max;   // buffer j end
        val[i][s_y - 1 - b][k] = rbuff_s[l] + v_min;   // buffer j start
      }
    }
  }
  
  /*----------------+
  |  K - direction  |
  +----------------*/
  if(dir == -1 || dir == 2) {

    /* not decomposed */
    if( dom->dim(Comp::k()) == 1 ) {
      if( bc().type(Dir::kmin(), BndType::periodic()) && 
          bc().type(Dir::kmax(), BndType::periodic()) ) {
        int kb_min = bc().index(Dir::kmin(), BndType::periodic());
        int kb_max = bc().index(Dir::kmax(), BndType::periodic());
        real v_min = bc().value(kb_min);
        real v_max = bc().value(kb_max);
        for_aij(i,j) {
          val[i][j][e_z + 1 + b] = val[i][j][s_z + o_z + b] + v_max;
          val[i][j][s_z - 1 - b] = val[i][j][e_z - o_z - b] + v_min;
        }
      }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_aij(i,j) {
        int l = j*ni()+i;
        assert(l < n*n);
        sbuff_e[l] = val[i][j][e_z - o_z - b];   // buffer k end
        sbuff_s[l] = val[i][j][s_z + o_z + b];   // buffer k start
        rbuff_e[l] = val[i][j][e_z +  1  + b];
        rbuff_s[l] = val[i][j][s_z -  1  - b];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], ni()*nj(),
                            par_real, dom->neighbour(Dir::kmax()),
                                      dom->neighbour(Dir::kmin()), Tag(4));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], ni()*nj(),
                            par_real, dom->neighbour(Dir::kmin()),
                                      dom->neighbour(Dir::kmax()), Tag(5));
#else
      for_aij(i,j) {
        int l = j*ni()+i;
        assert(l < n*n);
        rbuff_e[l] = val[i][j][e_z +  1 + b];
        rbuff_s[l] = val[i][j][s_z -  1 - b];
      }

      if( dom->neighbour(Dir::kmin()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_s[0], ni()*nj(), par_real,
                           dom->neighbour(Dir::kmin()), Tag(4), & req_r1 );
      }
      if( dom->neighbour(Dir::kmax()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_e[0], ni()*nj(), par_real,
                           dom->neighbour(Dir::kmax()), Tag(5), & req_r2 );
      }

      for_aij(i,j) {
        int l = j*ni()+i;
        assert(l < n*n);
        sbuff_e[l] = val[i][j][e_z - o_z - b];   // buffer k end
        sbuff_s[l] = val[i][j][s_z + o_z + b];   // buffer k start
      }

      if( dom->neighbour(Dir::kmax()) != par_proc_null ) {
        boil::cart.isend( &sbuff_e[0], ni()*nj(), par_real, 
                           dom->neighbour(Dir::kmax()), Tag(4), & req_s1 );
      }  
      if( dom->neighbour(Dir::kmin()) != par_proc_null ) {
        boil::cart.isend( &sbuff_s[0], ni()*nj(), par_real, 
                           dom->neighbour(Dir::kmin()), Tag(5), & req_s2 );
      }

      if( dom->neighbour(Dir::kmax()) != par_proc_null )
        boil::cart.wait( & req_s1, & ps );
      if( dom->neighbour(Dir::kmin()) != par_proc_null )
        boil::cart.wait( & req_s2, & ps );
      if( dom->neighbour(Dir::kmin()) != par_proc_null )
        boil::cart.wait( & req_r1, & ps );
      if( dom->neighbour(Dir::kmax()) != par_proc_null )
        boil::cart.wait( & req_r2, & ps );
#endif
      real v_min=0.0;
      real v_max=0.0;
      if(bc().type_here(Dir::kmin(), BndType::periodic())){
        int kb_min = bc().index(Dir::kmin(), BndType::periodic());
        v_min = bc().value(kb_min);
      }
      if(bc().type_here(Dir::kmax(), BndType::periodic())){
        int kb_max = bc().index(Dir::kmax(), BndType::periodic());
        v_max = bc().value(kb_max);
      }
      for_aij(i,j) {
        int l = j*ni()+i;
        assert(l < n*n);
        val[i][j][e_z + 1 + b] = rbuff_e[l] + v_max;   // buffer k end
        val[i][j][s_z - 1 - b] = rbuff_s[l] + v_min;   // buffer k start
      }
    }
  }

  }

  delete [] sbuff_s;
  delete [] sbuff_e;
  delete [] rbuff_s;
  delete [] rbuff_e;
}
