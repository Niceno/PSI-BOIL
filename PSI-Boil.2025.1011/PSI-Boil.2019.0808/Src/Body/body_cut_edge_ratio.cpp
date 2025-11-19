#include "body.h"

/******************************************************************************/
void Body::cut_edge_ratio(CutCell * ccell,Plane * p,
           const int n_in_fluid[][2][2],
           real x[], real y[], real z[]) {
 
    real er;

    /* cuts in "i" direction */
    for(int j=0; j<2; j++) {
      for(int k=0; k<2; k++) {
        er=-1.0;
        real xi = p->x( y[j], z[k] );
        if(n_in_fluid[0][j][k]==1 && n_in_fluid[1][j][k]==1){
          er=1.0;
        } else if(n_in_fluid[0][j][k]==0 && n_in_fluid[1][j][k]==0){
          er=0.0;
        } else if(n_in_fluid[0][j][k]==1 && n_in_fluid[1][j][k]==0){
          er = (xi-x[0])/(x[1]-x[0]);
        } else if(n_in_fluid[0][j][k]==0 && n_in_fluid[1][j][k]==1){
          er = (x[1]-xi)/(x[1]-x[0]);
        } else {
          if(n_in_fluid[0][j][k]==2){
            if(n_in_fluid[1][j][k]==0){
              er=0.0;
            } else if (n_in_fluid[1][j][k]==1){
              er=1.0;
            } else if (n_in_fluid[1][j][k]==2){
              er=0.0;
            }
          }
          if(n_in_fluid[1][j][k]==2){
            if(n_in_fluid[0][j][k]==0){
              er=0.0;
            } else if (n_in_fluid[0][j][k]==1){
              er=1.0;
            }
          }
        }
        assert( er!=-1.0 );
        ccell->fE(er,0,j,k);
      }
    }

    /* cuts in "j" direction */
    for(int i=0; i<2; i++){
      for(int k=0; k<2; k++) {
        er=-1.0;
        real yj = p->y( x[i], z[k] );
        if(n_in_fluid[i][0][k]==1 && n_in_fluid[i][1][k]==1){
          er=1.0;
        } else if(n_in_fluid[i][0][k]==0 && n_in_fluid[i][1][k]==0){
          er=0.0;
        } else if(n_in_fluid[i][0][k]==1 && n_in_fluid[i][1][k]==0){
            er = (yj-y[0])/(y[1]-y[0]);
        } else if(n_in_fluid[i][0][k]==0 && n_in_fluid[i][1][k]==1){
            er = (y[1]-yj)/(y[1]-y[0]);
        } else {
          if(n_in_fluid[i][0][k]==2){
            if(n_in_fluid[i][1][k]==0){
              er=0.0;
            } else if(n_in_fluid[i][1][k]==1){
              er=1.0;
            } else if(n_in_fluid[i][1][k]==2){
              er=0.0;
            }
          }
          if(n_in_fluid[i][1][k]==2){
            if(n_in_fluid[i][0][k]==0){
              er=0.0;
            } else if(n_in_fluid[i][0][k]==1){
              er=1.0;
            }
          }
        }
        assert( er!=-1.0 );
        ccell->fE(er,1,i,k);
      }
    }

    /* cuts in "k" direction */
    for(int i=0; i<2; i++) {
      for(int j=0; j<2; j++) {
        er=-1.0;
        real zk = p->z( x[i], y[j] );
        if (n_in_fluid[i][j][0]==1 && n_in_fluid[i][j][1]==1){
          er=1.0;
        } else if(n_in_fluid[i][j][0]==0 && n_in_fluid[i][j][1]==0){
          er=0.0;
        } else if(n_in_fluid[i][j][0]==1 && n_in_fluid[i][j][1]==0){
            er = (zk-z[0])/(z[1]-z[0]);
        } else if(n_in_fluid[i][j][0]==0 && n_in_fluid[i][j][1]==1){
            er = (z[1]-zk)/(z[1]-z[0]);
        } else {
          if(n_in_fluid[i][j][0]==2){
            if(n_in_fluid[i][j][1]==0){
              er=0.0;
            } else if(n_in_fluid[i][j][1]==1){
              er=1.0;
            } else if(n_in_fluid[i][j][1]==2){
              er=0.0;
            }
          }
          if(n_in_fluid[i][j][1]==2){
            if(n_in_fluid[i][j][0]==0){
              er=0.0;
            } else if(n_in_fluid[i][j][0]==1){
              er=1.0;
            }
          }
        }
        assert( er!=-1.0 );
        ccell->fE(er,2,i,j);
      }
    }
}
