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

  const real lx = grid_x_original->xn( grid_x_original->nnode() ) -
                  grid_x_original->xn( 1 );
  const real ly = grid_y_original->xn( grid_y_original->nnode() ) -
                  grid_y_original->xn( 1 );
  const real lz = grid_z_original->xn( grid_z_original->nnode() ) -
                  grid_z_original->xn( 1 );


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

  /* debugging info
  boil::aout << boil::pid << "ni =   " 
             << grid_x_local->ncell_b() << boil::endl;
  boil::aout << boil::pid << "nj =   " 
             << grid_y_local->ncell_b() << boil::endl;
  boil::aout << boil::pid << "nk =   " 
             << grid_z_local->ncell_b() << boil::endl;
  */

  // aout << "Decomposed the domain " << boil::endl;
  // aout << " c1x = " << cr_x.first() << boil::endl;
  // aout << " cNx = " << cr_x.last () << boil::endl;
  // aout << " c1y = " << cr_y.first() << boil::endl;
  // aout << " cNy = " << cr_y.last () << boil::endl;
  // aout << " c1z = " << cr_z.first() << boil::endl;
  // aout << " cNz = " << cr_z.last () << boil::endl;

  crsr = coarsen();

  //debug: if(boil::plot) boil::plot->plot(*this, name.c_str(), 0);
  
  boil::oout << "Domain decomposed into " << dim(Comp::i()) << " x " 
                                          << dim(Comp::j()) << " x "
                                          << dim(Comp::k()) << " blocks."
                                          << boil::endl;
}
