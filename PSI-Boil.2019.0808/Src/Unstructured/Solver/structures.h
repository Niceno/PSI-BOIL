#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "../../Global/global_constants.h"
#include "../../Global/global_minmax.h"

/*------+
|  ele  |
+------*/
struct ele {

  int  i,  j,  k, l;

  int  mark;
  int  n;          /* number of nodes */

  real x, y, Area; /* x, y, Area of a Delaunay cell */
 };

/*------+
|  sid  |
+------*/
struct sid {

  int  a, b;       /* left and right Delaunay cell */
  int  c, d;       /* start and end Voronoi cell */

  int  _c, _d;

  int  mark;       /* is on the boundary */

  real h, l;
 };

/*------+
|  nod  |
+------*/
struct nod {

  real x, y, Area; /* x, y, Area of Voronoi polygon) */

  int  mark;
  int  n;          /* number of neighbours */

  int  Ne;         /* number of elements */
  std::vector<int>  elem;
  std::vector<real> fe;
};

/*----------------------+
|  boundary conditions  |
+----------------------*/
struct bc {

  int  type;
  real u, v, T;

};

/*---------+
|  matrix  |
+---------*/
class matrix {

  public:
    int size() const {return val.size();}

    std::vector< std::vector<real> > val;
    std::vector< std::vector<int> >  con;
};

#endif
