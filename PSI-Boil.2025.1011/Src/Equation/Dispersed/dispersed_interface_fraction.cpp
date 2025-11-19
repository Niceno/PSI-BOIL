#include "dispersed.h"

/******************************************************************************/
real Dispersed::interface_fraction(int i, int j, int k, int p, 
                                   const real & cur_val) const {
/*-----------------------------------------------------------------------------+
|  if you set current cell value as parameter here, particles will not         |
|  overwrite each other.  If you send the value of "continuous", they will     |
|  both are needed!!!                                                          |
+-----------------------------------------------------------------------------*/
  real result = cur_val;

  const real xb = particles[p].x();
  const real yb = particles[p].y();
  const real zb = particles[p].z(); 
  const real radius = 0.5 * particles[p].d();

  /* use of sinusoidal interface */
  #if SIN
  const real dist = sqrt((xc(i)-xb)*(xc(i)-xb) + 
                         (yc(j)-yb)*(yc(j)-yb) + 
                         (zc(k)-zb)*(zc(k)-zb)) - radius; 

  /* this is a bit case-dependent, means a uniform mesh should be used */
  const real eps = 1.5 * dxc(boil::BW);

  if(dist < -eps) { 

     result = dispersed;

  } else if(dist < eps) { 

    real f_bub = (0.5 * (1.0 + (dist / eps) + 
                 (1.0 / boil::pi) * sin(boil::pi * dist / eps)));  

    const real f_dro = 1.0 - f_bub;
 
    /* don't erase the shape of current value (existing particle)  */
    /* volume fraction is defined always as the volume fraction of */ 
    /* the liquid in the cell                                      */
    if( dispersed < continuous ) { 
      result = boil::minr( f_bub, cur_val );
    } else if( dispersed > continuous ) { 
      result = boil::maxr( f_dro, cur_val );
    }   
  }

  /* use of volume fraction interface */
  #elif VF
  const real dist = sqrt((xc(i)-xb)*(xc(i)-xb) +
                         (yc(j)-yb)*(yc(j)-yb) +
                         (zc(k)-zb)*(zc(k)-zb)); 

  if (dist < (radius*0.75)) {

    result = dispersed;

  } else if(dist < (radius*1.25)) {

    int mm = 8;
    real  x0 = xn(i); real y0 = yn(j); real z0 = zn(k);
    real ddx = dxc(i)/real(mm);
    real ddy = dyc(j)/real(mm);
    real ddz = dzc(k)/real(mm);

    int itmp = 0;
    for (int ii = 0; ii < mm; ii++) {
      for (int jj = 0; jj < mm; jj++) {
        for (int kk = 0; kk < mm; kk++) {
          real xxc = x0 + 0.5 * ddx + real(ii)*ddx;
          real yyc = y0 + 0.5 * ddy + real(jj)*ddy;
          real zzc = z0 + 0.5 * ddz + real(kk)*ddz;
          real dist = pow(xxc-xb,2.0) + pow(yyc-yb,2.0) + pow(zzc-zb,2.0);
          if (dist > pow(radius,2.0)) itmp = itmp + 1;
        }   
      }   
    }   
    const real f_bub = real(itmp) /real(mm*mm*mm);
    const real f_dro = 1.0 - f_bub;

    /* don't erase the shape of current value (existing particle) */
    if(dispersed < continuous) {
      result = boil::minr( f_bub, cur_val );
    } else if( dispersed > continuous ) {
      result = boil::maxr( f_dro, cur_val );
    }
  }
  
  /* this is used in CIT */
  #elif CIT

    result = val[i][j][k];

  #endif

  return result;
}
