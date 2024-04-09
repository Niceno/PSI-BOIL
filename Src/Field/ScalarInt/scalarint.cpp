#include "scalarint.h"

/******************************************************************************/
ScalarInt::ScalarInt(const Domain & d) : alias(false) {
/*--------------------------------------+
|  creates a scalar for a given domain  |
+--------------------------------------*/

  dom = & d;
	
  allocate(d.ni(), d.nj(), d.nk());

  bndcnd = new BndCnd( *dom );

  nam = "";
}	

/******************************************************************************/
ScalarInt::ScalarInt(const Domain & d, const char * nm) : alias(false) {
/*--------------------------------------+
|  creates a scalar for a given domain  |
+--------------------------------------*/

  dom = & d;
	
  allocate(d.ni(), d.nj(), d.nk());

  bndcnd = new BndCnd( *dom );

  nam = nm;
}	

/******************************************************************************/
ScalarInt::ScalarInt(const Domain & d, BndCnd & b) : alias(false) {
/*--------------------------------------------------------------+
|  creates a scalar for a given domain and boundary conditions  |
+--------------------------------------------------------------*/

  dom = & d;

  allocate(d.ni(), d.nj(), d.nk());

  bndcnd = & b; 

  nam = "";
}	

/******************************************************************************/
ScalarInt::ScalarInt(const ScalarInt & s) : alias(false) {
/*---------------------------------------------------------------------+
|  creates a scalar like an existing scalar - used from Matrix object  |
+---------------------------------------------------------------------*/

  dom = s.domain();
	
  allocate(s.ni(), s.nj(), s.nk());

  bndcnd = new BndCnd( *dom ); 

  nam = "";
}	

/******************************************************************************/
ScalarInt::ScalarInt(const ScalarInt * s) : alias(true), nam(s->nam) {
/*------------------------------------+
|  this constructor creates an alias  |
+------------------------------------*/

  val    = s->val;

  n_x    = s->ni();
  s_x    = s->si();
  e_x    = s->ei();
  o_x    = s->ox();
  n_y    = s->nj();
  s_y    = s->sj();
  e_y    = s->ej();
  o_y    = s->oy();
  n_z    = s->nk();
  s_z    = s->sk();
  e_z    = s->ek();
  o_z    = s->oz();

  dom    = s->domain();
  bndcnd = s->bndcnd;
}	

/******************************************************************************/
void ScalarInt::allocate(int ni, int nj, int nk) {
 
  n_x = ni;
  n_y = nj;
  n_z = nk;
  s_x = boil::BW;
  s_y = boil::BW;
  s_z = boil::BW;
  e_x = ni-boil::BW-1;
  e_y = nj-boil::BW-1;
  e_z = nk-boil::BW-1;
  o_x = 0;
  o_y = 0;
  o_z = 0;

  alloc3d( &val, ni, nj, nk );
}

/******************************************************************************/
void ScalarInt::deallocate() {
/*-------------------------------------+
|  this is called from Vector as well  |
+-------------------------------------*/

  /* deallocate memory */
  dealloc3d( &val );
}	

/******************************************************************************/
ScalarInt::~ScalarInt() {

  if( !alias ) deallocate();
}
