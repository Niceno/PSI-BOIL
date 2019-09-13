#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
Domain::Domain(const Grid1D & ogx, const Grid1D & ogy, const Grid1D & ogz,
               const std::string n, const Decompose dec,
               const bool print_statistics) :
/*---------------------------------------------------------------+
|  this constructor creates a domain globally on all processors  |
+---------------------------------------------------------------*/
  lev(0), 
  cr_x(0,0), cr_y(0,0), cr_z(0,0),
  grid_x_original(&ogx), grid_y_original(&ogy), grid_z_original(&ogz),
  name(n) {

  body = new Empty();

  setup(dec);

  if(print_statistics)
    statistics();
}

/******************************************************************************/
Domain::Domain(const Grid1D & ogx, const Grid1D & ogy, const Grid1D & ogz,
               Body * stl_body,  
               const std::string n, const Decompose dec,
               const bool print_statistics) :
/*---------------------------------------------------------------+
|  this constructor creates a domain globally on all processors  |
+---------------------------------------------------------------*/
  lev(0), 
  cr_x(0,0), cr_y(0,0), cr_z(0,0),
  grid_x_original(&ogx), grid_y_original(&ogy), grid_z_original(&ogz),
  name(n) {

  body = stl_body;

  setup(dec);

  body->cut(*this);

  if(print_statistics)
    statistics(body);
}

/******************************************************************************/
Domain::Domain(const Grid1D * ogx, const Grid1D * ogy, const Grid1D * ogz,
               const Grid1D * lgx, const Grid1D * lgy, const Grid1D * lgz,
               int * dms, int * crds, int * nghbrs,
               const int l,  
               const Range<int> & crx, 
               const Range<int> & cry, 
               const Range<int> & crz,
               const std::string & n) :
  lev(l+1), cr_x(crx), cr_y(cry), cr_z(crz),
  grid_x_original(ogx), grid_y_original(ogy), grid_z_original(ogz), 
  grid_x_local(lgx),    grid_y_local(lgy),    grid_z_local(lgz),
  dims(dms), coords(crds), neighbours(nghbrs),
  name(n) {
/*-----------------------------------------------------------+
|  this constructor creates a domain locally on a processor  |
+-----------------------------------------------------------*/

  body = new Empty();

  crsr = coarsen();

  boil::oout << "Domain level " << level() << " created !" << boil::endl;
  //TMP if(boil::plot) boil::plot->plot(*this, name.c_str(), level());
}	
