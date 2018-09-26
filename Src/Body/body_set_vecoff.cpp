#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Field/ScalarBool/scalarbool.h"
#include "../Field/VectorBool/vectorbool.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::set_vecoff() {

  /* crude code!!! */
  /* vecoff should be set in accordance with geometry */

  for_m(m){
    int ii,jj,kk;
    ii=jj=kk=0;
    if(m==Comp::u()){
      ii=1;
    } else if (m==Comp::v()){
      jj=1;
    } else {
      kk=1;
    }

    for_vmijk((*vecoff), m, i, j, k){
      if( 0.5*(fV(i,j,k)+fV(i-ii,j-jj,k-kk)) <= 0.5 + boil::pico){
        (*vecoff)[m][i][j][k]=true;
      } else {
        (*vecoff)[m][i][j][k]=false;
      }
    }
    (*vecoff).exchange_all(m);

#if 0
    for_vmijk((*vecoff), m, i, j, k){
      if((*vecoff)[m][i][j][k]==true){
        fV(m,i,j,k)=0.0;
      }
    }
#endif
  }

}

/*-----------------------------------------------------------------------------+
 '$Id: body_set_vecoff.cpp,v 1.1 2014/02/03 14:12:34 sato Exp $'/
+-----------------------------------------------------------------------------*/
