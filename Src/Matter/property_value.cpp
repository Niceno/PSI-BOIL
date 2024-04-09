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
real Property::value(const Comp & m, 
                     const int i, const int j, const int k) const {

  if(con) return cval;
  int ox(0), oy(0), oz(0);
  if       (m==Comp::i()) {
    ox--;
  } else if(m==Comp::j()) {
    oy--;
  } else {
    oz--;
  }

  return 0.5*(val[i][j][k]+val[i+ox][j+oy][k+oz]);
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

/*============================================================================*/
real PropertyMix::value(const Comp & m,
                        const int i, const int j, const int k) const {
  int ox(0), oy(0), oz(0);
  if       (m==Comp::i()) {
    ox--;
  } else if(m==Comp::j()) {
    oy--;
  } else {
    oz--;
  }

  if(!bndc_a) {
    return 0.5*(value(i,j,k)+value(i+ox,j+oy,k+oz));
  } else {

    real col_a = (*bndc_a)[m][i][j][k];

    if( col_a > 1.0 ) col_a = 1.0;
    if( col_a < 0.0 ) col_a = 0.0;

    /* volume fraction of dispersed is treated using interpolation... */
    if( c_disp_a != NULL ) {
      real cdispa = 0.5*((*c_disp_a)[i][j][k]+(*c_disp_a)[i+ox][j+oy][k+oz]);
      if( cdispa > 0.0 ) 
         col_a = col_a + cdispa;
    }
    if( c_disp_b != NULL ) {
      real cdispb = 0.5*((*c_disp_b)[i][j][k]+(*c_disp_b)[i+ox][j+oy][k+oz]);
      if( cdispb < 1.0 ) 
           col_a = col_a + cdispb - 1.0;
    }

     return a -> value(m,i,j,k) * col_a +
            b -> value(m,i,j,k) * (1.0 - col_a);
  }
}
