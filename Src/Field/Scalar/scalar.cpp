#include "scalar.h"

/******************************************************************************/
Scalar::Scalar(const Domain & d) : alias(false) {
/*--------------------------------------+
|  creates a scalar for a given domain  |
+--------------------------------------*/

  dom = & d;
	
  allocate(d.ni(), d.nj(), d.nk());

  bndcnd = new BndCnd( *dom );

  nam = "";
}	

/******************************************************************************/
Scalar::Scalar(const Domain & d, const char * nm) : alias(false) {
/*--------------------------------------+
|  creates a scalar for a given domain  |
+--------------------------------------*/

  dom = & d;
	
  allocate(d.ni(), d.nj(), d.nk());

  bndcnd = new BndCnd( *dom );

  nam = nm;
}	

/******************************************************************************/
Scalar::Scalar(const Domain & d, BndCnd & b) : alias(false) {
/*--------------------------------------------------------------+
|  creates a scalar for a given domain and boundary conditions  |
+--------------------------------------------------------------*/

  dom = & d;

  allocate(d.ni(), d.nj(), d.nk());

  bndcnd = & b; 

  nam = "";
}	

/******************************************************************************/
Scalar::Scalar(const Scalar & s) : alias(false) {
/*---------------------------------------------------------------------+
|  creates a scalar like an existing scalar - used from Matrix object  |
+---------------------------------------------------------------------*/

  dom = s.domain();
	
  allocate(s.ni(), s.nj(), s.nk());

  bndcnd = new BndCnd( *dom ); 

  nam = "";
}	

/******************************************************************************/
Scalar::Scalar(const Scalar * s) : alias(true), nam(s->nam) {
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
void Scalar::allocate(int ni, int nj, int nk) {

  n_x = ni;   
  n_y = nj;   
  n_z = nk;
  s_x = 1;    
  s_y = 1;    
  s_z = 1;  
  e_x = ni-2; 
  e_y = nj-2; 
  e_z = nk-2; 
  o_x = 0;    
  o_y = 0;    
  o_z = 0;

  alloc3d( &val, ni, nj, nk );
}

/******************************************************************************/
void Scalar::deallocate() {
/*-------------------------------------+
|  this is called from Vector as well  |
+-------------------------------------*/

  /* deallocate memory */	
  dealloc3d( &val );
}	

/******************************************************************************/
Scalar::~Scalar() {

  if( !alias ) deallocate();
}	

/*-----------------------------------------------------------------------------+
 '$Id: scalar.cpp,v 1.25 2011/07/07 08:18:58 niceno Exp $'/
+-----------------------------------------------------------------------------*/
