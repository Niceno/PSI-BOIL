#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Body::cut_degen(real xp[], real yp[], real zp[], int & n,
                     int ic, int jc, int kc) {

  /* check degeneracy */
  int nda, ndb;
  int nnd;
  do {
    nnd=0;
    for(int nd1=0; nd1<=n-2; nd1++){
      for(int nd2=nd1+1; nd2<=n-1; nd2++){
        real dst = (xp[nd1]-xp[nd2])*(xp[nd1]-xp[nd2])
                  + (yp[nd1]-yp[nd2])*(yp[nd1]-yp[nd2])
                  + (zp[nd1]-zp[nd2])*(zp[nd1]-zp[nd2]);
        dst = sqrt(dst);
        if(dst<tol){
          nda=nd1;
          ndb=nd2;
          nnd++;
        }
      }
    }

    if(nnd>=1){
      for(int ntmp=ndb; ntmp<n; ntmp++){
        xp[ntmp]=xp[ntmp+1];
        yp[ntmp]=yp[ntmp+1];
        zp[ntmp]=zp[ntmp+1];
      }
      n--;
    }

  } while (nnd!=0);

}

/*-----------------------------------------------------------------------------+
 '$Id: body_cut_degen.cpp,v 1.1 2011/03/28 08:12:04 sato Exp $'/
+-----------------------------------------------------------------------------*/
