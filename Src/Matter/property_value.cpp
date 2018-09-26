#include "property.h"

/*============================================================================*/
real Property::value() const {

  assert(con==true);

  return cval;
}

/*============================================================================*/
real Property::value(const int i, const int j, const int k) const {
  if(con) return cval;
  return val[i][j][k];
}

/*============================================================================*/
void Property::value(const real & v) {

  assert(con==true);

  cval = v;
}

/*============================================================================*/
real Property::value_comp(const int comp) const {

  assert(comp >= 0 && comp <=1);

  if(comp == 1) return a -> value();
  else          return b -> value();
}

/*============================================================================*/
real PropertyMix::value(const int i, const int j, const int k) const {

  real col_a = (*c_a)[i][j][k];

  if( col_a > 1.0 ) col_a = 1.0;
  if( col_a < 0.0 ) col_a = 0.0;

  if( c_disp_a != NULL )
    if( (*c_disp_a)[i][j][k] > 0.0 ) 
       col_a = col_a + (*c_disp_a)[i][j][k];
//      col_a = boil::maxr( (*c_disp_a)[i][j][k], col_a ); 
  if( c_disp_b != NULL )
    if( (*c_disp_b)[i][j][k] < 1.0 ) 
         col_a = col_a + (*c_disp_b)[i][j][k] - 1.0;
//      col_a = boil::minr( (*c_disp_b)[i][j][k], col_a ); 

   return a -> value(i,j,k) * col_a +
          b -> value(i,j,k) * (1.0 - col_a);
}

/*-----------------------------------------------------------------------------+
 '$Id: property_value.cpp,v 1.3 2015/08/19 12:01:43 badreddine Exp $'/
+-----------------------------------------------------------------------------*/
