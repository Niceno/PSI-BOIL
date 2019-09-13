#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Domain::setup(const Decompose & dec) {
/*--------------------------------------------------------------+
|  this is called after global constructors, on all processors  |
+--------------------------------------------------------------*/

  if( (grid_x_original->periodic1() != grid_x_original->periodicN()) ||
      (grid_y_original->periodic1() != grid_y_original->periodicN()) || 
      (grid_z_original->periodic1() != grid_z_original->periodicN()) ) {
    boil::oout << "Inconsistent periodicity at Domain constructor" 
               << boil::endl;
    exit(0);
  }

  per[0] = grid_x_original->periodic1();
  per[1] = grid_y_original->periodic1();
  per[2] = grid_z_original->periodic1();

  ctf[0][0] = ( per[0]==0 && grid_x_original->bndgrid1() == BndGrid::symmetry() );
  ctf[0][1] = ( per[0]==0 && grid_x_original->bndgridN() == BndGrid::symmetry() );
  ctf[1][0] = ( per[1]==0 && grid_y_original->bndgrid1() == BndGrid::symmetry() );
  ctf[1][1] = ( per[1]==0 && grid_y_original->bndgridN() == BndGrid::symmetry() );
  ctf[2][0] = ( per[2]==0 && grid_z_original->bndgrid1() == BndGrid::symmetry() );
  ctf[2][1] = ( per[2]==0 && grid_z_original->bndgridN() == BndGrid::symmetry() );

  dummy[0] = grid_x_original->is_dummy();
  dummy[1] = grid_y_original->is_dummy();
  dummy[2] = grid_z_original->is_dummy();

#if 0
  const real lx = grid_x_original->xn( boil::BW + grid_x_original->ncell() ) -
                  grid_x_original->xn( boil::BW );
  const real ly = grid_y_original->xn( boil::BW + grid_y_original->ncell() ) -
                  grid_y_original->xn( boil::BW );
  const real lz = grid_z_original->xn( boil::BW + grid_z_original->ncell() ) -
                  grid_z_original->xn( boil::BW );
#else
  const real lx = grid_x_original->lx();
  const real ly = grid_y_original->lx();
  const real lz = grid_z_original->lx();
#endif

  /* initializes date for domain decomposition */
  /* (dims, coords, neighbours) */
  init(dec);
  
  /* get the local resolution ... */
  decompose(0, gi(), & cr_x);
  decompose(1, gj(), & cr_y);
  decompose(2, gk(), & cr_z);

  /* ... and create local grids */
  grid_x_local = new Grid1D( *grid_x_original, cr_x );
  grid_y_local = new Grid1D( *grid_y_original, cr_y );
  grid_z_local = new Grid1D( *grid_z_original, cr_z );

#if 0
  /* debugging info */
  boil::aout << boil::pid << "ni =   " 
             << grid_x_local->ncell_b() << boil::endl;
  boil::aout << boil::pid << "nj =   " 
             << grid_y_local->ncell_b() << boil::endl;
  boil::aout << boil::pid << "nk =   " 
             << grid_z_local->ncell_b() << boil::endl;
#endif

#if 0
  boil::aout << "Decomposed the domain " << boil::endl;
  boil::aout << " c1x = " << cr_x.first() << boil::endl;
  boil::aout << " cNx = " << cr_x.last () << boil::endl;
  boil::aout << " c1y = " << cr_y.first() << boil::endl;
  boil::aout << " cNy = " << cr_y.last () << boil::endl;
  boil::aout << " c1z = " << cr_z.first() << boil::endl;
  boil::aout << " cNz = " << cr_z.last () << boil::endl;
#endif

  crsr = coarsen();

  //debug: if(boil::plot) boil::plot->plot(*this, name.c_str(), 0);
  
  boil::oout << "Domain decomposed into " << dim(Comp::i()) << " x " 
                                          << dim(Comp::j()) << " x "
                                          << dim(Comp::k()) << " blocks."
                                          << boil::endl;
}
