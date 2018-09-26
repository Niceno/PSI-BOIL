#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "definitions.h"
#include "variables.h"

/******************************************************************************/
void form_Aphi() {

  for_sides(s) {
    const int c=side[s].c;
    const int d=side[s].d;

    Aphi.val[c][side[s]._d] = -side[s].l / side[s].h;
    Aphi.val[d][side[s]._c] = -side[s].l / side[s].h;

    Aphi.val[c][0] -= Aphi.val[c][side[s]._d];
    Aphi.val[d][0] -= Aphi.val[d][side[s]._c];
  }

  /*
  Aphi.val[Nn/2][0] *= 2.000;
  */
}

/*-----------------------------------------------------------------------------+
 '$Id: form_aphi.cpp,v 1.1 2015/09/01 09:41:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
