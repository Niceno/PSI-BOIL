#ifndef RANDOMFLOW_H
#define RANDOMFLOW_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include <vector>

#include "../Global/global_random.h"
#include "../Global/global_constants.h"
#include "../Parallel/Out/out.h"
#include "../Parallel/Out/print.h"
#include "../Field/Scalar/scalar.h"
#include "../Ravioli/comp.h"

////////////////////////////////////////////////
//                                            //
//  RandomFlow                                //
//  J Aerosol Science, Vol.25 (1994) pp91-112 //
//  Li et al., Aerosol ...                    //
//                                            //
////////////////////////////////////////////////
class RandomFlow {

  public:
    RandomFlow(int ne, real turb_time, real turb_length);
    void get_vector(const real t, const real *x,  real *v);
    void get_vector(const real t, const real *x,  real *v,
                    const Comp & i);

  private:
    real sclp(const real * A, const real * B); 
    void vecp(real * A, const real * B, const real * C); 
    real gauss();
    void gaussn(real **Y, real d, int n, int m); 
    void gaussn(real *Y, real d, int n); 

    int ne; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

    real
      *Omega,    /* Eq.15,celik.bib:\cite{LiAhetalJAS94} */
      **U1,**U2, /* velocity vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
      **K;       /* wave vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */

};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: randomflow.h,v 1.11 2015/01/05 16:28:18 sato Exp $'/
+-----------------------------------------------------------------------------*/
