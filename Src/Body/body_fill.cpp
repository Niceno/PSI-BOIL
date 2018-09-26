#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Body::push(int i, int j, int k) {
/*-------------------------+
|  just push on the stack  |
+-------------------------*/

  assert(stack_pointer < stack_size);

  stack[stack_pointer][0] = i;
  stack[stack_pointer][1] = j;
  stack[stack_pointer][2] = k;
  stack_pointer++;
}    

/******************************************************************************/
bool Body::pop(int * i, int * j, int * k) {
/*----------------------------------------------------------+
|  get from the stack if available, otherwise return false  |
+----------------------------------------------------------*/

  if(stack_pointer > 0) {
    stack_pointer--;
    * i = stack[stack_pointer][0];
    * j = stack[stack_pointer][1];
    * k = stack[stack_pointer][2];
    return true;
  }    
  else
   return false;
}   
 
/******************************************************************************/
void Body::body_fills(real new_color, int i, int j, int k) {
/*-----------------------------------------+
|  this is a three-dimensional flood-fill  |
|   (it puts one wherever it finds zero)   |
+-----------------------------------------*/
 
  do { /* do ... */

    (*sca)[i][j][k] = new_color;

    if(i+1 <= eifl[3] && (*sca)[i+1][j][k] < -0.5) push(i+1,j,k);
    if(i-1 >= sifl[3] && (*sca)[i-1][j][k] < -0.5) push(i-1,j,k);

    if(j+1 <= ejfl[3] && (*sca)[i][j+1][k] < -0.5) push(i,j+1,k);
    if(j-1 >= sjfl[3] && (*sca)[i][j-1][k] < -0.5) push(i,j-1,k);

    if(k+1 <= ekfl[3] && (*sca)[i][j][k+1] < -0.5) push(i,j,k+1);
    if(k-1 >= skfl[3] && (*sca)[i][j][k-1] < -0.5) push(i,j,k-1);

  }
  while( pop(&i,&j,&k) ); /* ... while there are cells on the stack */
}

/******************************************************************************/
void Body::body_fillv(real new_color, int i, int j, int k,const Comp m) {
/*-----------------------------------------+
|  this is a three-dimensional flood-fill  |
|   (it puts one wherever it finds zero)   |
+-----------------------------------------*/

  do { /* do ... */

    (*vec)[m][i][j][k] = new_color;

    if(i+1 <= (*vec).ei(m) && (*vec)[m][i+1][j][k] < -0.5) push(i+1,j,k);
    if(i-1 >= (*vec).si(m) && (*vec)[m][i-1][j][k] < -0.5) push(i-1,j,k);

    if(j+1 <= (*vec).ej(m) && (*vec)[m][i][j+1][k] < -0.5) push(i,j+1,k);
    if(j-1 >= (*vec).sj(m) && (*vec)[m][i][j-1][k] < -0.5) push(i,j-1,k);

    if(k+1 <= (*vec).ek(m) && (*vec)[m][i][j][k+1] < -0.5) push(i,j,k+1);
    if(k-1 >= (*vec).sk(m) && (*vec)[m][i][j][k-1] < -0.5) push(i,j,k-1);

  }
  while( pop(&i,&j,&k) ); /* ... while there are cells on the stack */
}

/*-----------------------------------------------------------------------------+
 '$Id: body_fill.cpp,v 1.7 2011/03/28 08:12:04 sato Exp $'/
+-----------------------------------------------------------------------------*/
