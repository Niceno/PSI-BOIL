#include "scalar.h"

/******************************************************************************/
void Scalar::one2all(const int sender) {
	
  if( boil::cart.nproc() < 2 ) return;

  assert( sender < boil::cart.nproc() );

  const int size = ni() * nj() * nk();

  for(int receiver=0; receiver<boil::cart.nproc(); receiver++)
    if(receiver != sender) {

      /* check if the sizes match */
      int sizes[2] = {0, 0};
      if     (boil::cart.iam() == sender)   sizes[0] = size; /* sender size */
      else if(boil::cart.iam() == receiver) sizes[1] = size; /* reciver size */
      boil::cart.sum_int_n(sizes, 2);
      if(sizes[0] != sizes[1]) {
        APR(sizes[0]);
        APR(sizes[1]);
        boil::oout << "sizes do not match in one2all! exiting!" << boil::endl;
        exit(0);
      }

      /* if they do match, send the data */
      boil::cart.one2one(val[0][0], val[0][0], size, par_real, 
                         sender, receiver, Tag(7));
    }
}
