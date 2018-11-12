#include "scalarint.h"

/******************************************************************************/
void ScalarInt::bnd_extract( const Dir d, int *** cplane ) const { 

  //std::cout<<"scalar_bnd_extract:begin\n";
  /* buffer for parallel runs and temporary plane */
  int *  buffer;
  int ** tplane;

  /* number of subdomains in j and k directions */
  const int di = dom->dim(Comp::i());
  const int dj = dom->dim(Comp::j());
  const int dk = dom->dim(Comp::k());

  /*-----------------------------------+
  |  direction is either imin or imax  |
  +-----------------------------------*/
  if( (d == Dir::imin()) || (d == Dir::imax()) ) {

    /* allocate memory for the copy plane */
    alloc2d(   cplane, nj() , nk() );

    /* browses through quadrants of domain decomposition */
    for(int dcj = 0; dcj < dj; dcj++) 
      for(int dck = 0; dck < dk; dck++) { 

        /* allocate memory for the buffers (this also zeroes them) */
        alloc2d( &tplane, nj() , nk() );
        alloc1d( &buffer, nj() * nk() );

        if( (dom->coord(Comp::j()) == dcj) && (dom->coord(Comp::k()) == dck) ) {
          if( (d == Dir::imin()) && (dom->coord(Comp::i()) == 0) ) 
            for_ajk(j, k) 
              tplane[j][k] = val[si()-1][j][k];
              //tplane[j][k] = val[si()][j][k];

          if( (d == Dir::imax()) && (dom->coord(Comp::i()) == di-1) )
            for_ajk(j, k) 
              tplane[j][k] = val[ei()+1][j][k];
              //tplane[j][k] = val[ei()][j][k];
        } 

        /* copy the plane into the buffer (note that, at this point, all
           domains but the one close to the boundaries will have zeroes */
        int c = 0;
        for_ajk(j, k) 
          buffer[c++] = tplane[j][k];

        /* distribute the buffer over all processors */
        boil::cart.sum_int_n(buffer, nj()*nk());

        /* take the value only if in the proper quadrant */
        if( (dom->coord(Comp::j()) == dcj) && (dom->coord(Comp::k()) == dck) ) {
          int c = 0;
          for_ajk(j, k) 
            (*cplane)[j][k] = buffer[c++];
        } 

        dealloc1d( &buffer );
        dealloc2d( &tplane );
      }
  }

  /*-----------------------------------+
  |  direction is either jmin or jmax  |
  +-----------------------------------*/
  if( (d == Dir::jmin()) || (d == Dir::jmax()) ) {

    /* allocate memory for the copy plane and the buffer */
    alloc2d(  cplane, ni() , nk() );

    for(int dci = 0; dci < di; dci++) 
      for(int dck = 0; dck < dk; dck++) { 

        /* allocate memory for the buffers (this also zeroes them) */
        alloc2d( &tplane, ni() , nk() );
        alloc1d( &buffer, ni() * nk() );

        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::k()) == dck) ) {
          if( (d == Dir::jmin()) && (dom->coord(Comp::j()) == 0) ) 
            for_aik(i, k) 
              //tplane[i][k] = val[i][sj()-1][k];
              tplane[i][k] = val[i][sj()][k];

          if( (d == Dir::jmax()) && (dom->coord(Comp::j()) == dj-1) )
            for_aik(i, k) 
              //tplane[i][k] = val[i][ej()+1][k];
              tplane[i][k] = val[i][ej()][k];
        }

        /* copy the plane into the buffer (note that, at this point, all
           domains but the one close to the boundaries will have zeroes */
        int c = 0;
        for_aik(i, k) 
          buffer[c++] = tplane[i][k];

        /* distribute the buffer over all processors */
        boil::cart.sum_int_n(buffer, ni()*nk());

        /* take the value only if in the proper quadrant */
        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::k()) == dck) ) {
          int c = 0;
          for_aik(i, k) 
            (*cplane)[i][k] = buffer[c++];
        } 

        dealloc1d( &buffer );
        dealloc2d( &tplane );
      }
  }

  /*-----------------------------------+
  |  direction is either kmin or kmax  |
  +-----------------------------------*/
  if( (d == Dir::kmin()) || (d == Dir::kmax()) ) {

    /* allocate memory for the copy plane and the buffer */
    alloc2d(  cplane, ni() , nj() );

    for(int dci = 0; dci < di; dci++) 
      for(int dcj = 0; dcj < dj; dcj++) { 

        /* allocate memory for the buffers (this also zeroes them) */
        alloc2d( &tplane, ni() , nj() );
        alloc1d( &buffer, ni() * nj() );

        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::j()) == dcj) ) {
          if( (d == Dir::kmin()) && (dom->coord(Comp::k()) == 0) ) 
            for_aik(i, j) 
              //tplane[i][j] = val[i][j][sk()-1];
              tplane[i][j] = val[i][j][sk()];

          if( (d == Dir::kmax()) && (dom->coord(Comp::k()) == dk-1) )
            for_aik(i, j) 
              //tplane[i][j] = val[i][j][ek()+1];
              tplane[i][j] = val[i][j][ek()];
        }

        /* copy the plane into the buffer (note that, at this point, all
           domains but the one close to the boundaries will have zeroes */
        int c = 0;
        for_aij(i, j) 
          buffer[c++] = tplane[i][j];

        /* distribute the buffer over all processors */
        boil::cart.sum_int_n(buffer, ni()*nj());

        /* take the value only if in the proper quadrant */
        if( (dom->coord(Comp::i()) == dci) && (dom->coord(Comp::j()) == dcj) ) {
          int c = 0;
          for_aij(i, j) 
            (*cplane)[i][j] = buffer[c++];
        } 

        dealloc1d( &buffer );
        dealloc2d( &tplane );
      }   
  }

  //std::cout<<"scalar_bnd_extract:end\n";
}
