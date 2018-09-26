#include "body.h"

/******************************************************************************/
void Polygon::arrange() {

   /*-----------------------------------------+ 
   |  take absolute values of surface normal  |
   +-----------------------------------------*/
   const real anx = fabs(nor[0]);
   const real any = fabs(nor[1]);
   const real anz = fabs(nor[2]);

   const real maxan = boil::maxr( anx, any, anz );

   /*--------------------+
   |  center of gravity  |
   +--------------------*/
   real g[3] = {0.0, 0.0, 0.0}; /* gx, gy, gz */
   for(int n=0; n<nn; n++) {
     g[0] += x[n];
     g[1] += y[n];
     g[2] += z[n];
   }
   for(int i=0; i<3; i++) g[i] /= (real)nn;
     
   std::vector<AngleNode> nodes;

   /*-----------------------+
   |  sort nodes by angles  |
   +-----------------------*/
     
   real alfa;
   /* place all nodes into array "nodes" */
   for(int n=0; n<nn; n++) {
     if( maxan == anz ) alfa = nor[2]/anz * atan2( (y[n]-g[1]), (x[n]-g[0]) );
     if( maxan == any ) alfa = nor[1]/any * atan2( (x[n]-g[0]), (z[n]-g[2]) );
     if( maxan == anx ) alfa = nor[0]/anx * atan2( (z[n]-g[2]), (y[n]-g[1]) );
            
     AngleNode an( alfa, x[n], y[n], z[n]);
     nodes.push_back(an);
   }

   /* sort the array "nodes" */
   stable_sort(nodes.begin(), nodes.end());

   /* place the sorted vales from "nodes" back to coordinates */
   for(int n=0; n<nn; n++) {
     x[n] = nodes[n].x;
     y[n] = nodes[n].y;
     z[n] = nodes[n].z;
   }
}

/*-----------------------------------------------------------------------------+
 '$Id: polygon_arrange.cpp,v 1.1 2015/08/17 10:37:55 niceno Exp $'/
+-----------------------------------------------------------------------------*/
