#include "grid1d.h"

/******************************************************************************/
Grid1D::Grid1D(const Range<real> &  xr, const Range<real> & dxr,
               const int & n, const Periodic & p) : N(n), 
  period1(p), 
  periodN(p) {
/*---------------------------+
|  creates non-uniform grid  |
+---------------------------*/

  allocate(); 

  distribute_nodes_inside( xr.first(),  xr.last(), 
                          dxr.first(), dxr.last(), N);
  
  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Range<real> &  xr, 
               const int & n, const Periodic & p) : N(n), 
               period1(p), periodN(p) {
/*-----------------------+
|  creates uniform grid  | -> derived from above; could only one be used???
+-----------------------*/

  allocate(); 

  real dx = (xr.last() - xr.first()) / (real)N;

  distribute_nodes_inside( xr.first(),  xr.last(), dx, dx, N);
  
  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D & left, const Grid1D & right, 
               const Periodic & p) : N(left.ncell()+right.ncell()), 
  period1(p), 
  periodN(p) {
/*---------------------------+
|  creates non-uniform grid  |
+---------------------------*/

  allocate(); 

  /* just copty the first grid */
  for(int i=1; i<=left.nnode(); i++)
    x_node[i] = left.xn(i);

  /* append the second grid */ 
  for(int i=2; i<=right.nnode(); i++)
    x_node[i+left.ncell()] = x_node[left.nnode()] + right.xn(i) - right.xn(1);

  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D & left, const Grid1D & center, const Grid1D & right, 
               const Periodic & p) : 
  N(left.ncell()+center.ncell()+right.ncell()), 
  period1(p), 
  periodN(p) {
/*---------------------------+
|  creates non-uniform grid  |
+---------------------------*/

  allocate(); 

  /* just copty the first grid */
  for(int i=1; i<=left.nnode(); i++)
    x_node[i] = left.xn(i);

  /* append the second grid */ 
  for(int i=2; i<=center.nnode(); i++)
    x_node[i+left.ncell()] = x_node[left.nnode()] + center.xn(i) - center.xn(1);

  /* append the third grid */ 
  for(int i=2; i<=right.nnode(); i++)
    x_node[i+left.ncell() + center.ncell()] 
    = 
    x_node[left.nnode() + center.nnode() - 1] + right.xn(i) - right.xn(1);

  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D & grid, 
               const Step   & step) // default step is Step(1)
  : period1(grid.periodic1()), 
    periodN(grid.periodicN()) {
/*-------------------------------------------------------------------------+
|  copy constructor with possibility to coarsen. not sure if it is needed  | 
|                                                                          |
|                        N = 8                                             |
|                                                                          |
|        0   1   2   3   4   5   6   7   8   9                             |
|      | o |-O-+-O-+-O-+-O-+-O-+-O-+-O-+-O-| o |                           |
|      0   1   2   3   4   5   6   7   8   9  10                           |
|                                                                          |
|                        N = 4                                             |
|                                                                          |
|      0       1       2       3       4       5                           |
|  |   o   |---O---+---O---+---O---+---O---|   o   |                       |
|  0       1       2       3       4       5       6                       |
+-------------------------------------------------------------------------*/

  /* check if it is possible */
  if(grid.ncell() % step.size() !=0) {
    boil::aout << "Critical! Can not coarsen grid! Stopping!" << boil::endl;
    exit(0);
  }

  /* if yes, create (coarser) grid */
  N = grid.ncell()/step.size();

  allocate(); 

  for(int i=1; i<N+2; i++)
    x_node[i] = grid.xn(i*step.size()-step.size()+1);

  int i=N+2;
  x_node[i] = 2.0 * x_node[i-1] - x_node[i-2];
  
  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D     & grid, 
               const Range<int> & cr,     // first and last cell inside
               const Step       & step) 
  : period1(grid.periodic1()), 
    periodN(grid.periodicN()) {
/*-------------------------------------------------------------------------+
|  copy constructor which gets a subgrid.                                  | 
|  the dissatvantage is that it allocates and de-allocates memory :-(      |
+-------------------------------------------------------------------------*/

  /* check if it is possible */
  if( (cr.last() - cr.first() + 1) % step.size() != 0 ) {
    boil::aout << "Critical! Can not coarsen grid! Stopping!" << boil::endl;
    exit(0);
  }

  /* make a temporary coarser grid */
  Grid1D * coarse = new Grid1D(grid, step);

  /* handle periodicity */
  if(cr.first() > 1)            period1 = Periodic::no(); 
  if(cr.last()  < grid.ncell()) periodN = Periodic::no(); 

  /* create cell1 and cellN for coarse grid */
  const int cell1 = (cr.first()+step.size()-1)/step.size();
  const int cellN =  cr.last()/step.size();

  /* set the right number of cells */
  N = cellN - cell1 + 1;

  allocate(); 

  /* just copy the values */
  for(int i=0; i<nnode_b(); i++) {
    const int j = i+cell1-1;
     x_node[i] = coarse-> xn(j);
    dx_node[i] = coarse->dxn(j);
  }
  for(int i=0; i<ncell_b(); i++) {
    const int j = i+cell1-1;
     x_cell[i] = coarse-> xc(j);
    dx_cell[i] = coarse->dxc(j);
  }

  /* clean up */
  delete coarse;
}  

/******************************************************************************/
void Grid1D::allocate() {

   x_cell = new real[ncell_b()];
  dx_cell = new real[ncell_b()];
   x_node = new real[nnode_b()];
  dx_node = new real[nnode_b()];
}

/******************************************************************************/
Grid1D::~Grid1D() {
  N = 0;
  delete []  x_node;
  delete []  x_cell;
  delete [] dx_cell;
  delete [] dx_node;
}	

/******************************************************************************/
void Grid1D::print() const {

  boil::aout << boil::endl;
  boil::aout << "-[ node based dimensions ]-------------------------------" 
             << boil::endl;
  if(periodic1()==true) boil::aout << ":::::::::::" << boil::endl;
  for(int i=0; i<nnode_b(); i++) {
    if(i == 1 || i == nnode_b()-1) boil::aout << "- - - - - -" << boil::endl; 
    boil::aout << "i="<<i << "   xn="<<xn(i) << "   dxn="<<dxn(i) 
               << boil::endl;
  }

  if(periodicN()==true) boil::aout << ":::::::::::" << boil::endl;
  boil::aout << boil::endl;
  boil::aout << "-[ cell based dimensions ]-------------------------------" 
             << boil::endl;
  if(periodic1()==true) boil::aout << ":::::::::::" << boil::endl;
  for(int i=0; i<ncell_b(); i++) {
    if(i == 1 || i == ncell_b()-1) boil::aout << "- - - - - -" << boil::endl; 
    boil::aout << "i="<<i << "   xc="<<xc(i) << "   dxc="<<dxc(i) 
               << boil::endl;
  }
  if(periodicN()==true) boil::aout << ":::::::::::" << boil::endl;
}

/******************************************************************************/
void Grid1D::plot(const char * name) const {

  Board board;

  /* grid dimensions and thick width */
  const real l  = xn(nnode_b()-2) - xn(1);
  const real hw = l/100.0;
  const real rc = hw*0.5;

  board.setLineWidth( 1.0 );

  /*---------------------+
  |  plot inside domain  |
  +---------------------*/
  board.setPenColorRGBi(  0, 0,  0);

  /* domain */
  board.drawLine( xn(0), 0, xn(nnode_b()-1), 0 );

  /* nodes */
  for(int i=0; i<ncell_b(); i++) {
    board.drawLine( xn(i), -hw, xn(i), hw );
  }

  /* cells */
  for(int i=1; i<ncell_b()-1; i++) {
    
    const real xc = 0.5*(xn(i)+xn(i+1));
    board.fillCircle( xc, 0, rc );
  }

  /*----------------+
  |  plot boundary  |
  +----------------*/
  board.setPenColorRGBi(255,  0,  0);

  /* domain */
  board.drawLine( xn(0),            0 , xn(1),           0  );
  board.drawLine( xn(nnode_b()-1),  0 , xn(nnode_b()-2), 0  );

  /* nodes */
  board.drawLine( xn(0),           -hw, xn(0),           hw );
  board.drawLine( xn(nnode_b()-1), -hw, xn(nnode_b()-1), hw );

  /* cells */
  board.fillCircle( 0.5*(xn(0)          +xn(1)          ), 0, rc );
  board.fillCircle( 0.5*(xn(nnode_b()-1)+xn(nnode_b()-2)), 0, rc );

  board.saveEPS( name );
}
/*-----------------------------------------------------------------------------+
 '$Id: grid1d.cpp,v 1.21 2016/02/10 08:38:23 sato Exp $'/
+-----------------------------------------------------------------------------*/
