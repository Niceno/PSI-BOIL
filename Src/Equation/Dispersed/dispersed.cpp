#include "dispersed.h"
#include "dispersed_browsing.h"

/******************************************************************************
  This is a non-trivial constructor. 
*******************************************************************************/
Dispersed::Dispersed(const Scalar & c, 
                           Scalar * color, 
                           int disp, 
                     const Vector & v,
                     const Times  & t, 
                     const Matter * f,
                     const Matter * s) :
 
  Scalar(& c), /* it creates an alias, not a new field!!! */
  col   (color), 
  p_id (*(c.domain())), 
  u     (& v),
  time (& t),
  flu  (f),
  sol  (s),
  dom(c.domain()),

  box_diam_ratio(1.5),  /* used for particle box */ 

  uvw_limit(0.35),

  list_diameter(0.004), /* used in cell_list.cpp as grid cell size: for collisions */

  continuous(1-disp), 

  dispersed(disp) {                    
 
  p_id = color->shape();
  assert(col == color);
  
  cell_init();

  *this = continuous; /* sets "color function" to continuous */
  p_id = OFF + boil::milli;
 
}

/******************************************************************************/
Dispersed::~Dispersed() {
/*--------------------------------------------------------------------------+
|  destructor is also placed here.                                          |
|  this is a psi-boil-style comment, written with lower-case letters only.  |
+--------------------------------------------------------------------------*/

}

/*-----------------------------------------------------------------------------+
 '$Id: dispersed.cpp,v 1.9 2015/08/19 11:39:11 badreddine Exp $'/ 
+-----------------------------------------------------------------------------*/
