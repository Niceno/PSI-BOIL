#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Domain::decompose(const int i, const int g, Range<int> * cr) const {

  /* array which will store number of cells in this direction */
  /* it can be used later to find the offset numbers for each direction */
  int * n_dir = new int [ dims[i] ];

  /* number of global cells inside (without halo/border cells) */
  const int g_in = g-2*boil::BW; 

  /* minimum number of cells in each processor */
  const int n_min = g_in / dims[i];

  for(int k=0; k<dims[i]; k++) 
    n_dir[k] = n_min;

  /* remaining cells */
  int n_res = g_in - n_min * dims[i];
  if(n_res != 0) 
    boil::oout << "Warning: domains will not be balanced!!!" << boil::endl;

  for(int k=dims[i]-1; k>=dims[i]-n_res; k--) 
    n_dir[k]++;

  /* debugging info
  for(int k=0; k<dims[i]; k++)
    boil::aout << boil::pid << "n_dir[k] = " << n_dir[k] << boil::endl;
  */  

  /* estimate cell range */
  int n_left_side = 0;
  for(int k=0; k<coords[i]; k++)
    n_left_side += n_dir[ k ];

  cr->first( n_left_side + 1 );
  cr->last ( n_left_side + n_dir[ coords[i] ] );

  delete [] n_dir;
}
