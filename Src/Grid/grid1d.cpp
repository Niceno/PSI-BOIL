#include "grid1d.h"

/******************************************************************************/
Grid1D::Grid1D(const Range<real> &  xr, const Range<real> & dxr,
               const int & n, const Periodic & p) : nc_in(n), 
  period1(p), 
  periodN(p) {
/*---------------------------+
|  creates non-uniform grid  |
+---------------------------*/

  allocate(); 

  distribute_nodes_inside( xr.first(),  xr.last(),
                          dxr.first(), dxr.last());
  
  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Range<real> &  xr,
               const int & n, const Periodic & p) : nc_in(n),  
               period1(p), periodN(p) {
/*-----------------------+
|  creates uniform grid  | -> derived from above; could only one be used???
+-----------------------*/

  allocate(); 

  real dx = (xr.last() - xr.first()) / (real)nc_in;

  distribute_nodes_inside( xr.first(),  xr.last(), dx, dx);

  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D & left, const Grid1D & right,
               const Periodic & p) : nc_in(left.ncell()+right.ncell()), 
               period1(p), periodN(p) {
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
  nc_in(left.ncell()+center.ncell()+right.ncell()), 
  period1(p), periodN(p) {
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
|                    nc_in = 8                                             |
|                                                                          |
|        0   1   2   3   4   5   6   7   8   9                             |
|      | o |-O-+-O-+-O-+-O-+-O-+-O-+-O-+-O-| o |                           |
|      0   1   2   3   4   5   6   7   8   9  10                           |
|                                                                          |
|                    nc_in = 4                                             |
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
  nc_in = grid.ncell()/step.size();

  allocate(); 

  for(int i=0; i<=nc_in; i++)
    x_node[boil::BW + i] = grid.xn(boil::BW + i*step.size());

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
  nc_in = cellN - cell1 + 1;

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
  nc_in = 0;
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

