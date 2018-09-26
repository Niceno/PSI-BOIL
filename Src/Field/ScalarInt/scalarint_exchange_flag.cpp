#include "scalarint.h"
//#define DEBUG
/* uncomment the line below for blocking (old) send-receive */
//#define SENDRECV

/******************************************************************************/
void ScalarInt::exchange(const int * ical, const int dir) const {

#ifdef DEBUG
  std::cout<<"exchange_flag::start "<<boil::cart.iam()<<" "
           <<ical[boil::cart.iam()]<<"\n";
#endif

  if( bc().count() == 0 ) {
    OMS(Warning: exchanging a variable without boundary conditions);
  }

  /*-------------------------------+
  |                                |
  |  buffers for parallel version  |
  |                                |
  +-------------------------------*/
  int * sbuff_s, * sbuff_e, * rbuff_s, * rbuff_e;

  par_request req_s1,req_s2,req_r1,req_r2;
  par_status  ps;

  /* allocate memory for buffers */
  const int n = boil::maxi(ni(), nj(), nk())+1;
  assert(n > 0);

  sbuff_s = new int [n*n];
  sbuff_e = new int [n*n];
  rbuff_s = new int [n*n];
  rbuff_e = new int [n*n];

  int irank = boil::cart.iam();

  /*----------------+
  |  I - direction  |
  +----------------*/
  if(dir == -1 || dir == 0) {

    /* not decomposed in I direction */
    if( dom->dim(Comp::i()) == 1 ) {
      if( bc().type(Dir::imin(), BndType::periodic()) && 
          bc().type(Dir::imax(), BndType::periodic()) )
        for_jk(j,k) {
          val[e_x+1][j][k] = val[s_x + o_x][j][k];
          val[s_x-1][j][k] = val[e_x - o_x][j][k];
        }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_jk(j,k) {
        int l = k*nj()+j;
        sbuff_e[l] = val[e_x - o_x][j][k];   // buffer i end
        sbuff_s[l] = val[s_x + o_x][j][k];   // buffer i start
        rbuff_e[l] = val[e_x +  1 ][j][k];
        rbuff_s[l] = val[s_x -  1 ][j][k];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], nj()*nk(),
                            par_int, dom->neighbour(Dir::imax()),
                                      dom->neighbour(Dir::imin()), Tag(0));
        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], nj()*nk(),
                            par_int, dom->neighbour(Dir::imin()),
                                      dom->neighbour(Dir::imax()), Tag(1));
#else
      for_jk(j,k) {
        int l = k*nj()+j;
        rbuff_e[l] = val[e_x +  1 ][j][k];
        rbuff_s[l] = val[s_x -  1 ][j][k];
      }

      if( dom->neighbour(Dir::imin()) != par_proc_null ) {
        if( ical[dom->neighbour(Dir::imin())] !=0 ){
          //std::cout<<"exchange_flag::imin receive "<<irank<<"\n";
          boil::cart.irecv( &rbuff_s[0], nj()*nk(), par_int, 
                           dom->neighbour(Dir::imin()), Tag(0), & req_r1 );
      } }
      if( dom->neighbour(Dir::imax()) != par_proc_null ) {
        if( ical[dom->neighbour(Dir::imax())] !=0 ){
          //std::cout<<"exchange_flag::imax receive "<<irank<<"\n";
          boil::cart.irecv( &rbuff_e[0], nj()*nk(), par_int, 
                           dom->neighbour(Dir::imax()), Tag(1), & req_r2 );
      } }

      for_jk(j,k) {
        int l = k*nj()+j;
        sbuff_e[l] = val[e_x - o_x][j][k];   // buffer i end
        sbuff_s[l] = val[s_x + o_x][j][k];   // buffer i start
      }

      if( dom->neighbour(Dir::imax()) != par_proc_null ) {
        if( ical[irank] !=0 ){
          //std::cout<<"exchange_flag::imax send "<<irank<<"\n";
          boil::cart.isend( &sbuff_e[0], nj()*nk(), par_int, 
                           dom->neighbour(Dir::imax()), Tag(0), & req_s1 );
      } }
      if( dom->neighbour(Dir::imin()) != par_proc_null ) {
        if( ical[irank] !=0 ){
          //std::cout<<"exchange_flag::imin send "<<irank<<"\n";
          boil::cart.isend( &sbuff_s[0], nj()*nk(), par_int, 
                           dom->neighbour(Dir::imin()), Tag(1), & req_s2 );
      } }

      if( dom->neighbour(Dir::imax()) != par_proc_null ){
        if( ical[irank] !=0 ){
          boil::cart.wait( & req_s1, & ps );
          //std::cout<<"exchange_flag::wait:imax send "<<irank<<"\n";
      } }
      if( dom->neighbour(Dir::imin()) != par_proc_null ){
        if( ical[irank] !=0 ){
          boil::cart.wait( & req_s2, & ps );
          //std::cout<<"exchange_flag::wait:imin send "<<irank<<"\n";
      } }
      if( dom->neighbour(Dir::imin()) != par_proc_null ){
        if( ical[dom->neighbour(Dir::imin())] !=0 ){
          boil::cart.wait( & req_r1, & ps );
          //std::cout<<"exchange_flag::wait:imin receive "<<irank<<"\n";
      } }
      if( dom->neighbour(Dir::imax()) != par_proc_null ){
        if( ical[dom->neighbour(Dir::imax())] !=0 ){
          boil::cart.wait( & req_r2, & ps );
          //std::cout<<"exchange_flag::wait:imax receive "<<irank<<"\n";
      } }
#endif
      for_jk(j,k) {
        int l = k*nj()+j;
        val[e_x+1][j][k] = rbuff_e[l];   // buffer i end
        val[s_x-1][j][k] = rbuff_s[l];   // buffer i start
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
          bc().type(Dir::jmax(), BndType::periodic()) )
        for_ik(i,k) {
          val[i][e_y+1][k] = val[i][s_y + o_y][k];
          val[i][s_y-1][k] = val[i][e_y - o_y][k];
        }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_ik(i,k) {
        int l = k*ni()+i;
        sbuff_e[l] = val[i][e_y - o_y][k];   // buffer j end
        sbuff_s[l] = val[i][s_y + o_y][k];   // buffer j start
        rbuff_e[l] = val[i][e_y +  1 ][k];
        rbuff_s[l] = val[i][s_y -  1 ][k];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], ni()*nk(),
                            par_int, dom->neighbour(Dir::jmax()),
                                      dom->neighbour(Dir::jmin()), Tag(2));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], ni()*nk(),
                            par_int, dom->neighbour(Dir::jmin()),
                                      dom->neighbour(Dir::jmax()), Tag(3));
#else
      for_ik(i,k) {
        int l = k*ni()+i;
        rbuff_e[l] = val[i][e_y +  1 ][k];
        rbuff_s[l] = val[i][s_y -  1 ][k];
      }

      if( dom->neighbour(Dir::jmin()) != par_proc_null ) {
        if( ical[dom->neighbour(Dir::jmin())] !=0 ){
          boil::cart.irecv( &rbuff_s[0], ni()*nk(), par_int, 
                             dom->neighbour(Dir::jmin()), Tag(2), & req_r1 );
      } }
      if( dom->neighbour(Dir::jmax()) != par_proc_null ) {
        if( ical[dom->neighbour(Dir::jmax())] !=0 ){
          boil::cart.irecv( &rbuff_e[0], ni()*nk(), par_int, 
                             dom->neighbour(Dir::jmax()), Tag(3), & req_r2 );
      } }

      for_ik(i,k) {
        int l = k*ni()+i;
        sbuff_e[l] = val[i][e_y - o_y][k];   // buffer j end
        sbuff_s[l] = val[i][s_y + o_y][k];   // buffer j start
      }
  
      if( dom->neighbour(Dir::jmax()) != par_proc_null ) {
        if( ical[irank] !=0 ){
          boil::cart.isend( &sbuff_e[0], ni()*nk(), par_int, 
                             dom->neighbour(Dir::jmax()), Tag(2), & req_s1 );
      } }
      if( dom->neighbour(Dir::jmin()) != par_proc_null ) {
        if( ical[irank] !=0 ){
          boil::cart.isend( &sbuff_s[0], ni()*nk(), par_int, 
                             dom->neighbour(Dir::jmin()), Tag(3), & req_s2 );
      } }

      if( dom->neighbour(Dir::jmax()) != par_proc_null )
        if( ical[irank] !=0 )
          boil::cart.wait( & req_s1, & ps );
      if( dom->neighbour(Dir::jmin()) != par_proc_null )
        if( ical[irank] !=0 )
          boil::cart.wait( & req_s2, & ps );
      if( dom->neighbour(Dir::jmin()) != par_proc_null )
        if( ical[dom->neighbour(Dir::jmin())] !=0 )
          boil::cart.wait( & req_r1, & ps );
      if( dom->neighbour(Dir::jmax()) != par_proc_null )
        if( ical[dom->neighbour(Dir::jmax())] !=0 )
          boil::cart.wait( & req_r2, & ps );
#endif
      for_ik(i,k) {
        int l = k*ni()+i;
        val[i][e_y+1][k] = rbuff_e[l];   // buffer j end
        val[i][s_y-1][k] = rbuff_s[l];   // buffer j start
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
          bc().type(Dir::kmax(), BndType::periodic()) )
        for_ij(i,j) {
          val[i][j][e_z+1] = val[i][j][s_z + o_z];
          val[i][j][s_z-1] = val[i][j][e_z - o_z];
        }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_ij(i,j) {
        int l = j*ni()+i;
        sbuff_e[l] = val[i][j][e_z - o_z];   // buffer k end
        sbuff_s[l] = val[i][j][s_z + o_z];   // buffer k start
        rbuff_e[l] = val[i][j][e_z +  1 ];
        rbuff_s[l] = val[i][j][s_z -  1 ];
      
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], ni()*nj(),
                            par_int, dom->neighbour(Dir::kmax()),
                                      dom->neighbour(Dir::kmin()), Tag(4));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], ni()*nj(),
                            par_int, dom->neighbour(Dir::kmin()),
                                      dom->neighbour(Dir::kmax()), Tag(5));
#else
      for_ij(i,j) {
        int l = j*ni()+i;
        rbuff_e[l] = val[i][j][e_z +  1 ];
        rbuff_s[l] = val[i][j][s_z -  1 ];
      }

      if( dom->neighbour(Dir::kmin()) != par_proc_null ) {
        if( ical[dom->neighbour(Dir::kmin())] !=0 ){
          boil::cart.irecv( &rbuff_s[0], ni()*nj(), par_int,
                             dom->neighbour(Dir::kmin()), Tag(4), & req_r1 );
      } }
      if( dom->neighbour(Dir::kmax()) != par_proc_null ) {
        if( ical[dom->neighbour(Dir::kmax())] !=0 ){
        boil::cart.irecv( &rbuff_e[0], ni()*nj(), par_int,
                           dom->neighbour(Dir::kmax()), Tag(5), & req_r2 );
      } }

      for_ij(i,j) {
        int l = j*ni()+i;
        sbuff_e[l] = val[i][j][e_z - o_z];   // buffer k end
        sbuff_s[l] = val[i][j][s_z + o_z];   // buffer k start
      }

      if( dom->neighbour(Dir::kmax()) != par_proc_null ) {
        if( ical[irank] !=0 ){
          boil::cart.isend( &sbuff_e[0], ni()*nj(), par_int, 
                             dom->neighbour(Dir::kmax()), Tag(4), & req_s1 );
      } } 
      if( dom->neighbour(Dir::kmin()) != par_proc_null ) {
        if( ical[irank] !=0 ){
        boil::cart.isend( &sbuff_s[0], ni()*nj(), par_int, 
                           dom->neighbour(Dir::kmin()), Tag(5), & req_s2 );
      } }

      if( dom->neighbour(Dir::kmax()) != par_proc_null )
        if( ical[irank] !=0 )
          boil::cart.wait( & req_s1, & ps );
      if( dom->neighbour(Dir::kmin()) != par_proc_null )
        if( ical[irank] !=0 )
          boil::cart.wait( & req_s2, & ps );
      if( dom->neighbour(Dir::kmin()) != par_proc_null )
        if( ical[dom->neighbour(Dir::kmin())] !=0 )
          boil::cart.wait( & req_r1, & ps );
      if( dom->neighbour(Dir::kmax()) != par_proc_null )
        if( ical[dom->neighbour(Dir::kmax())] !=0 )
          boil::cart.wait( & req_r2, & ps );
#endif
      for_ij(i,j) {
        int l = j*ni()+i;
        val[i][j][e_z+1] = rbuff_e[l];   // buffer k end
        val[i][j][s_z-1] = rbuff_s[l];   // buffer k start
      }
    }
  }
  
  delete [] sbuff_s;
  delete [] sbuff_e;
  delete [] rbuff_s;
  delete [] rbuff_e;
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarint_exchange_flag.cpp,v 1.1 2015/05/05 14:36:01 sato Exp $'/
+-----------------------------------------------------------------------------*/
