#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
const Domain * Domain::coarsen() const {
/*---------------------------------------------+ 
|  create grids for the next (coarser) level   |
+---------------------------------------------*/

  /* to disable coarsening, uncomment the following line */
  //return NULL;
 
  /* minimum resolution (with buffers) */
  const int min_n =  std::max(4,boil::BW) + 2*boil::BW;

  /* coarsening factors */
  int c_fac_x = 2; 
  int c_fac_y = 2; 
  int c_fac_z = 2; 

  /* accumulated steps (multiplication of coarsening factors) */
  /*------------------------------------------+ 
  |  be very careful with this since static   |
  |  variables are shared between objects of  |
  |             the same class!!!             |
  +------------------------------------------*/
  static int step_x = 1;
  static int step_y = 1;
  static int step_z = 1;

  /* estimate factors (only first one is important) */
  int c_factors[64], nf;
#if 0
  factor( ni()-2*boil::BW, c_factors, &nf); if(nf>0) c_fac_x = c_factors[0];
  factor( nj()-2*boil::BW, c_factors, &nf); if(nf>0) c_fac_y = c_factors[0];
  factor( nk()-2*boil::BW, c_factors, &nf); if(nf>0) c_fac_z = c_factors[0];
#else
  if(!grid_x_local->is_dummy()) {
    factor( ni()-2*boil::BW, c_factors, &nf); 
    if(nf>0) c_fac_x = c_factors[0];
  }
  if(!grid_y_local->is_dummy()) {
    factor( nj()-2*boil::BW, c_factors, &nf); 
    if(nf>0) c_fac_y = c_factors[0];
  }
  if(!grid_z_local->is_dummy()) {
    factor( nk()-2*boil::BW, c_factors, &nf);
    if(nf>0) c_fac_z = c_factors[0];
  }
#endif

  const Grid1D * g_x_coarse = grid_x_local;
  const Grid1D * g_y_coarse = grid_y_local;
  const Grid1D * g_z_coarse = grid_z_local;

  int ijk_crble[] = { 1, 1, 1 }; // coarsenable

  /* only domain */
  if(    ((ni() - 2*boil::BW) % c_fac_x!=0)
      ||  (ni() <= min_n) 
      ||   ni() - 2*boil::BW == c_fac_x
    ) ijk_crble[0] = 0;
  if(    ((nj() - 2*boil::BW) % c_fac_y!=0) 
      ||  (nj() <= min_n) 
      ||   nj() - 2*boil::BW == c_fac_y
    ) ijk_crble[1] = 0;
  if(   ((nk() - 2*boil::BW) % c_fac_z!=0) 
      || (nk() <= min_n) 
      ||  nk() - 2*boil::BW == c_fac_z
    ) ijk_crble[2] = 0;

  /* gather information over other domains */
  boil::cart.sum_int_n(ijk_crble, 3);

  if( ijk_crble[0] == boil::cart.nproc() ) {
    step_x *= c_fac_x;
    g_x_coarse = new Grid1D( *grid_x_original, cr_x, Step( step_x ));
  }

  if( ijk_crble[1] == boil::cart.nproc() ) {
    step_y *= c_fac_y;
    g_y_coarse = new Grid1D( *grid_y_original, cr_y, Step( step_y ));
  }

  if( ijk_crble[2] == boil::cart.nproc() ) {
    step_z *= c_fac_z;
    g_z_coarse = new Grid1D( *grid_z_original, cr_z, Step( step_z ));
  }


  /* if any direction is coarsened */
  if( (g_x_coarse != grid_x_local) || 
      (g_y_coarse != grid_y_local) ||
      (g_z_coarse != grid_z_local) ) {
    return new Domain( grid_x_original, grid_y_original, grid_z_original,   
                       g_x_coarse, g_y_coarse, g_z_coarse,
                       dims, coords, neighbours,
                       level(), cr_x, cr_y, cr_z,
                       name );
  } 
  /* else */

  /*--------------------------------------------+ 
  |  important: re-initialize static variables  |
  +--------------------------------------------*/
  step_x = 1;
  step_y = 1;
  step_z = 1;

  return NULL;
}	
