#include "dispersed.h"

/******************************************************************************/
void Dispersed::box_velocity() {

  /*==================================================================+
  |                                                                   |
  |  browse through all particles to store velocities in each corner  |
  |                                                                   |
  +==================================================================*/

  for_p(p) { /* browse through all particles of the "disp" class */

    /*--------------------------------------+
    |  browse through all corners of a box  |
    |                                       |
    |            6---------------7          |
    |           /|              /|          |
    |       ke 4---------------5 |          |
    |          | |             | |          |
    |          | |             | |          |
    |          | |             | |          |
    |          | 2-------------|-3 je       |
    |          |/              |/           |
    |       ks 0---------------1 js         |
    |         is              ie            |
    +--------------------------------------*/
    int corner = 0;
    for(int kk=0; kk<2; kk++) {
      for(int jj=0; jj<2; jj++) {
        for(int ii=0; ii<2; ii++) {

          int i = particles[p].si();
          int j = particles[p].sj();
          int k = particles[p].sk();
          if(ii>0) i = particles[p].ei();
          if(jj>0) j = particles[p].ej();
          if(kk>0) k = particles[p].ek();

          /* set velocities to zero */
          real uvwc[DIM] = {0.0, 0.0, 0.0}; 

          /* fetch the central value */ 
           if(dom->contains_ijk(i,j,k)) {
              u->central(i, j, k, uvwc);
           }

          /* store velocity. note it will be non-zero only 
             in processors which contain the I,J,K coordinate */
          particles[p].box_u(corner, uvwc[~Comp::u()]);
          particles[p].box_v(corner, uvwc[~Comp::v()]);
          particles[p].box_w(corner, uvwc[~Comp::w()]);

          /* increase the corner count (important) */
          corner++;
        }
      }
    }
  }

  /*=====================================+
  |                                      |
  |  exchange values between processors  |
  |                                      |
  +=====================================*/
  const int buffer_size = particles.size() * 8 * DIM;
  real * buffer     = new real[buffer_size];
  for(int ww = 0; ww < buffer_size; ww++) buffer[ww] = 0.0;

  /*-----------------------------------------+ 
  |  copy the values from box_uvw to buffer  |
  +-----------------------------------------*/ 
  int count = 0;
  for_dp((*this), p) {
    for(int corner = 0; corner < 8; corner++) {
      buffer[count++] = particles[p].box_u(corner);
      buffer[count++] = particles[p].box_v(corner);
      buffer[count++] = particles[p].box_w(corner);
    }
  }
  assert(count == buffer_size);
 
  /*------------------------------+ 
  |  summ up over all processors  |
  +------------------------------*/ 
  boil::cart.sum_real_n(buffer, buffer_size);

  /* at this point, all processors have all 
     velocities at all corners of all particle boxes */

  /*-----------------------------------------+ 
  |  copy the values from buffer to box_uvw  |
  +-----------------------------------------*/ 
  count = 0;
  for_dp((*this), p) 
    for(int corner = 0; corner < 8; corner++) {
      particles[p].box_u(corner, buffer[count++]);
      particles[p].box_v(corner, buffer[count++]);
      particles[p].box_w(corner, buffer[count++]);
    }
  assert(count == buffer_size);

  delete [] buffer;

}
