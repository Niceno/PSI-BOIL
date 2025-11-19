#include "body.h"

/******************************************************************************/
real distance(const real x0, const real x1, 
              const real y0, const real y1,
              const real z0, const real z1) {

  return sqrt(  ( x0 - x1 ) * ( x0 - x1 ) +
                ( y0 - y1 ) * ( y0 - y1 ) +
                ( z0 - z1 ) * ( z0 - z1 ) ) ;
}

/******************************************************************************/
void Polygon::split(Polygon ** poly_left,
                    Polygon ** poly_right) const {
 
  std::vector<int> nodes_left;      
  std::vector<int> nodes_right;

  real xl[6];
  real yl[6];
  real zl[6];

  real xr[6];
  real yr[6];
  real zr[6];

  const int no = nnodes();

  if(no == 3) return;

  /* create combinations */
  int comb[9][2] = { {0,2}, {0,3}, {0,4},
                     {1,3}, {1,4}, {1,5},
                     {2,4}, {2,5},       
                     {3,5} }; 
  const int n_comb = 9;

  /* initialize distances */
  real dist[9];
  for(int c=0; c<n_comb; c++) 
    dist[c] = boil::giga;
  

  for(int c=0; c<n_comb; c++) {
    const int na = comb[c][0];
    const int nb = comb[c][1];

    /* take onodes_lefty possible combinations */
    if( (nb < no) && (nb < (na+no-1))  ) {
      dist[c] = distance(x[na], x[nb], y[na], y[nb], z[na], z[nb]);
    }

    /* if hexagon, allow only big diagonals */
    if( no == 6 ) {
      if( (nb - na) != 3 )
        dist[c] = boil::giga;
    }
  }

  /* find the smallest distances */
  real sd = dist[0];
  int  id = 0; 
  for(int c=0; c<n_comb; c++) {
    if( dist[c] < sd ) {
      sd = dist[c];
      id = c; 
    }
  }

  /* fill up the left polygon */
  int i = comb[id][0];
  do {
    nodes_left.push_back(i++);
    if(i == no) 
      i = 0; 
  }
  while(i != comb[id][1]);
  nodes_left.push_back(i);

  /* fill up the right polygon */
  int j = comb[id][1];
  do {
    nodes_right.push_back(j++);
    if(j == no) 
      j = 0; 
  }
  while(j != comb[id][0]);
  nodes_right.push_back(j);

  /* handle/copy coordinates */
  for(int c=0; c<nodes_left.size(); c++) {
    xl[c] = x[ nodes_left[c] ];
    yl[c] = y[ nodes_left[c] ];
    zl[c] = z[ nodes_left[c] ];
  }
  for(int c=0; c<nodes_right.size(); c++) {
    xr[c] = x[ nodes_right[c] ];
    yr[c] = y[ nodes_right[c] ];
    zr[c] = z[ nodes_right[c] ];
  }

  /* create new polygons */
  * poly_left  = new Polygon( (int)nodes_left.size(),  xl, yl, zl );
  * poly_right = new Polygon( (int)nodes_right.size(), xr, yr, zr );
}
