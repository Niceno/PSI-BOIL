#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "variables.h"
#include "definitions.h"

/******************************************************************************/
void calc_con() {
/*--------------------------------+
|  find the connectivity matrices |
+--------------------------------*/

  /* For delaunay nodes (voronoi centers) */
  for_sides(s) {
    const int c=side[s].c;
    const int d=side[s].d;

    Auvw.con[c][0]++;
    Auvw.con[c][Auvw.con[c][0]]=d;
    side[s]._d = Auvw.con[c][0];

    Auvw.con[d][0]++;
    Auvw.con[d][Auvw.con[d][0]]=c;
    side[s]._c = Auvw.con[d][0];

    Aphi.con[c][0]++;
    Aphi.con[c][Aphi.con[c][0]]=d;
    side[s]._d = Aphi.con[c][0];

    Aphi.con[d][0]++;
    Aphi.con[d][Aphi.con[d][0]]=c;
    side[s]._c = Aphi.con[d][0];
  }
}
