#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Domain::distribute(const int nproc, const int dim, int * dis, int * res) {

  /*---------------------------+
  |  initialize distributions  |
  +---------------------------*/
  for(int m=0; m<dim; m++)
    dis[m] = 1;
  
  /*---------------------------------+
  |  factorize number of processors  |
  +---------------------------------*/
  int factors[64], nfact;
  factor( nproc, factors, &nfact);      	

  /*-------------------------------------------------+
  |  assign factor to coordinate directions (dis[])  |
  +-------------------------------------------------*/
  for(int f=0; f<nfact; f++) {
  
    int maxr = boil::maxi( res[0], res[1], res[2] );

    bool done = false;

    /* try the biggest one */
    for(int d=0; d<dim; d++) {
      if( maxr == res[d] ) {
        if( res[d]%factors[f] == 0 ) {
          dis[d] *= factors[f];
          res[d] /= factors[f];
          done = true;
          break;
        }
      }
    }

    /* if it fails, go for the one which is dividable */
    if( !done ) {
      for(int d=0; d<dim; d++) {
        if( res[d]%factors[f] == 0 ) {
          dis[d] *= factors[f];
          res[d] /= factors[f];
          done = true;
          break;
        }
      }
    }

    assert(done);
  }
}
 
/*-----------------------------------------------------------------------------+
 '$Id: domain_distribute.cpp,v 1.4 2009/07/17 16:36:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
