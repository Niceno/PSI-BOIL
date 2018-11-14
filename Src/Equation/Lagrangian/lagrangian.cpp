#include "lagrangian.h"
#include "lagrangian_browsing.h"

/******************************************************************************
  This is a non-trivial constructor. 
*******************************************************************************/
Lagrangian::Lagrangian(const Scalar & c, 
                             Scalar * color, 
                             real     lagp,
                       const Scalar & cfv,
                       const Vector & v,
                       const Times  & t, 
                       const Matter * f,
                       const Matter * s) :
 
  Scalar(& c),  // it creates an alias, not a new field!!!
  col   (color), 
  p_id (*(c.domain())), 
  cfu  (& cfv),
  u    (& v),
  time (& t),
  flu  (f),
  sol  (s),

  dom(c.domain()),

  box_diam_ratio(1.5),  /* used for particle box */ 

  //uvw_limit(0.35),  //not-used

  list_diameter(0.004), /* used in cell_list.cpp as grid cell size: for collisions */

  continuous(1-lagp),  //not-really-in-use

  lagrangian(lagp) {
    assert(lagrangian == 0.56789);  //mark

    p_id = color->shape();
    assert(col == color);
  
    cell_init();

    *this = continuous; /* sets "color function" to continuous */
    p_id = OFF + boil::milli;
 
}

/******************************************************************************/
Lagrangian::~Lagrangian() {
/*--------------------------------------------------------------------------+
|  destructor is also placed here.                                          |
|  this is a psi-boil-style comment, written with lower-case letters only.  |
+--------------------------------------------------------------------------*/

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian.cpp,v 1.9 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
