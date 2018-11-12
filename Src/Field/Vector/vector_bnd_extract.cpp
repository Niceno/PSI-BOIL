#include "vector.h"

/******************************************************************************/
void Vector::bnd_extract( const Dir d, real **** cplane ) const { 

  /* buffer for parallel runs and temporary plane */
  real *  buffer;
  real *** tplane;

  /* number of subdomains in j and k directions */
  const int di = dom->dim(Comp::i());
  const int dj = dom->dim(Comp::j());
  const int dk = dom->dim(Comp::k());

  /*-----------------------------------+
  |  direction is either imin or imax  |
  +-----------------------------------*/
  if( (d == Dir::imin()) || (d == Dir::imax()) ) {

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, nj()+1, nk()+1);

    /* browses through quadrants of domain decomposition */
    for(int dcj = 0; dcj < dj; dcj++)
      for(int dck = 0; dck < dk; dck++) {

        /* allocate memory for the buffers (this also zeroes them) */
        alloc3d( &tplane, 3, nj()+1 , nk()+1 );
        alloc1d( &buffer, 3 * (nj()+1) * (nk()+1) );

        if( (dom->coord(Comp::j()) == dcj) && (dom->coord(Comp::k()) == dck) ) {
          if( (d == Dir::imin()) && (dom->coord(Comp::i()) == 0) )
            for_m(m){
              //std::cout<<"bnd_extract:"<<m<<" "<<si(m)<<"\n";
              for_amjk(m, j, k) 
                tplane[~m][j][k] = vec[m][si(m)][j][k];
            }

          if( (d == Dir::imax()) && (dom->coord(Comp::i()) == di-1) )
            for_m(m){
              //std::cout<<"bnd_extract:"<<m<<" "<<ei(m)<<"\n";
              for_amjk(m, j, k) {
                tplane[~m][j][k] = vec[m][ei(m)][j][k];
              }
            }
          //exit(0);
        }

        /* copy the plane into the buffer (note that, at this point, all
           domains but the one close to the boundaries will have zeroes */
        int c = 0;
        for_m(m){
          for_amjk(m, j, k){
            buffer[c++] = tplane[~m][j][k];
          }
        }

        /* distribute the buffer over all processors */
        boil::cart.sum_real_n(buffer, 3*nj()*nk());

        /* take the value only if in the proper quadrant */
        if( (dom->coord(Comp::j()) == dcj) && (dom->coord(Comp::k()) == dck) ) {
          int c = 0;
          for_m(m)
            for_amjk(m, j, k)
              (*cplane)[~m][j][k] = buffer[c++];
        }

        dealloc1d( &buffer );
        dealloc3d( &tplane );
      }
  }

  /*-----------------------------------+
  |  direction is either jmin or jmax  |
  +-----------------------------------*/
  if( (d == Dir::jmin()) || (d == Dir::jmax()) ) {

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, ni()+1, nk()+1);

    /* browses through quadrants of domain decomposition */
    for(int dci = 0; dci < di; dci++)
      for(int dck = 0; dck < dk; dck++) {

        /* allocate memory for the buffers (this also zeroes them) */
        alloc3d( &tplane, 3, ni()+1 , nk()+1 );
        alloc1d( &buffer, 3 * (ni()+1) * (nk()+1) );

        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::k()) == dck) ) {
          if( (d == Dir::jmin()) && (dom->coord(Comp::j()) == 0) )
            for_m(m)
              for_amik(m, i, k) 
                tplane[~m][i][k] = vec[m][i][sj(m)][k];

          if( (d == Dir::jmax()) && (dom->coord(Comp::j()) == dj-1) )
            for_m(m)
              for_amik(m, i, k) 
                tplane[~m][i][k] = vec[m][i][ej(m)][k];
        }

        /* copy the plane into the buffer (note that, at this point, all
           domains but the one close to the boundaries will have zeroes */
        int c = 0;
        for_m(m)
          for_amik(m, i, k)
            buffer[c++] = tplane[~m][i][k];

        /* distribute the buffer over all processors */
        boil::cart.sum_real_n(buffer, 3*ni()*nk());

        /* take the value only if in the proper quadrant */
        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::k()) == dck) ) {
          int c = 0;
          for_m(m)
            for_amik(m, i, k)
              (*cplane)[~m][i][k] = buffer[c++];
        }

        dealloc1d( &buffer );
        dealloc3d( &tplane );
      }
  }

  /*-----------------------------------+
  |  direction is either kmin or kmax  |
  +-----------------------------------*/
  if( (d == Dir::kmin()) || (d == Dir::kmax()) ) {

    /* allocate memory for the copy plane */
    alloc3d( cplane, 3, ni()+1, nj()+1);

    /* browses through quadrants of domain decomposition */
    for(int dci = 0; dci < di; dci++)
      for(int dcj = 0; dcj < dj; dcj++) {

        /* allocate memory for the buffers (this also zeroes them) */
        alloc3d( &tplane, 3, ni()+1 , nj()+1 );
        alloc1d( &buffer, 3 * (ni()+1) * (nj()+1) );

        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::j()) == dcj) ) {
          if( (d == Dir::kmin()) && (dom->coord(Comp::k()) == 0) )
            for_m(m)
              for_amij(m, i, j) 
                tplane[~m][i][j] = vec[m][i][j][sk(m)];

          if( (d == Dir::kmax()) && (dom->coord(Comp::k()) == dk-1) )
            for_m(m)
              for_amij(m, i, j) 
                tplane[~m][i][j] = vec[m][i][j][ek(m)];
        }

        /* copy the plane into the buffer (note that, at this point, all
           domains but the one close to the boundaries will have zeroes */
        int c = 0;
        for_m(m)
          for_amij(m, i, j)
            buffer[c++] = tplane[~m][i][j];

        /* distribute the buffer over all processors */
        boil::cart.sum_real_n(buffer, 3*ni()*nj());

        /* take the value only if in the proper quadrant */
        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::j()) == dcj) ) {
          int c = 0;
          for_m(m)
            for_amij(m, i, j)
              (*cplane)[~m][i][j] = buffer[c++];
        }

        dealloc1d( &buffer );
        dealloc3d( &tplane );
      }
  }
}
