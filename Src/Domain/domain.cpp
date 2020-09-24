#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

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
  name(n), dc(dec) {

  if(stl_body) {
    body = stl_body;
  } else {
    body = new Empty();
  }

  setup(dc);

  if(stl_body) {
    body->cut(*this);

    if(print_statistics)
      statistics(body);
  } else {
    if(print_statistics)
      statistics();
  }
}

/******************************************************************************/
Domain::Domain(const Domain & fine_dom, 
               const Step cx, const Step cy, const Step cz,
               Body * stl_body,  
               const bool print_statistics) :
/*---------------------------------------------------------------+
|  this constructor creates a domain globally on all processors  |
+---------------------------------------------------------------*/
  lev(0),
  cr_x(0,0), cr_y(0,0), cr_z(0,0),
  name(fine_dom.dom_name()+" - coarsened"), 
  dc(fine_dom.decomp()) {

  const Step * c1 = &cx;  
  const Step * c2 = &cy;  
  const Step * c3 = &cz;  

  if       (cy.size()<0&&cz.size()<0) {
    c2 = &cx;
    c3 = &cx;
  } else if( cx.size()<=0 || cy.size()<=0 || cz.size()<=0 ) {
    boil::oout<<"Disallowed step values! "<<cx<<" "<<cy<<" "<<cz<<" .\n"
              <<"Exiting."<<boil::endl;
    exit(0);
  }

  if(!boil::cart.iam())
    boil::oout<<"Creating coarser grid with strides "
              <<*c1<<" "<<*c2<<" "<<*c3<<" .\n";

  /* properly coarsen dummy grids */
  if(fine_dom.grid_x_org()->is_dummy()) {
    grid_x_original = new Grid1D(c1->size()*fine_dom.grid_x_org()->lx());
  } else {
    grid_x_original = new Grid1D(*fine_dom.grid_x_org(),*c1);
  }

  if(fine_dom.grid_y_org()->is_dummy()) {
    grid_y_original = new Grid1D(c2->size()*fine_dom.grid_y_org()->lx());
  } else {
    grid_y_original = new Grid1D(*fine_dom.grid_y_org(),*c2);
  }

  if(fine_dom.grid_z_org()->is_dummy()) {
    grid_z_original = new Grid1D(c3->size()*fine_dom.grid_z_org()->lx());
  } else {
    grid_z_original = new Grid1D(*fine_dom.grid_z_org(),*c3);
  }

  if(stl_body) {
    body = stl_body;
  } else {
    body = new Empty();
  }

  setup(dc);

  if(stl_body) {
    body->cut(*this);

    if(print_statistics)
      statistics(body);
  } else {
    if(print_statistics)
      statistics();
  }
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
  name(n), dc(Decompose::undefined()) {
/*-----------------------------------------------------------+
|  this constructor creates a domain locally on a processor  |
+-----------------------------------------------------------*/

  body = new Empty();

  crsr = coarsen();

  boil::oout << "Domain level " << level() << " created !" << boil::endl;
  //TMP if(boil::plot) boil::plot->plot(*this, name.c_str(), level());
}	

/******************************************************************************/
/*----------------+
|  static members |
+----------------*/
int Domain::factor_x = 1;
int Domain::factor_y = 1;
int Domain::factor_z = 1;

void Domain::set_decomposition_factors(const int fx,
                                       const int fy,
                                       const int fz) {
  factor_x = fx;
  factor_y = fy;
  factor_z = fz;

  return;
} 
