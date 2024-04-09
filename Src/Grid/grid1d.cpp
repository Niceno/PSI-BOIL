#include "grid1d.h"

/* If the cutoff at N is undefined, the same is used for both sides */
/* BndGrids are only relevant for non-periodic grids */
/******************************************************************************/
Grid1D::Grid1D(const Range<real> &  xr, const Range<real> & dxr,
               const int & n, const Periodic & p, 
               const BndGrid & co1, const BndGrid & coN) : 
  nc_in(n),
  nc_tot(n+2*boil::BW),
  dummy_grid(false),
  period1(p), 
  periodN(p),
  ctf1(co1),
  ctfN(coN) {

#if 1
   if(nc_in<boil::BW) {
     boil::aout<<"At least as many cells as the buffer width ("<<boil::BW
               <<") are required in each direction. Exiting."<<boil::endl;
     exit(0);
   }
#endif

/*----------------+
|  check cutoffs  |
+----------------*/
  if(p == Periodic::no()&&ctf1 == BndGrid::undefined()) {
    boil::aout<<"At least one cutoff must be specified for a non-periodic grid."
              <<" Exiting."<<boil::endl;
    exit(0);
  }
  if(p == Periodic::no()&&ctfN == BndGrid::undefined()) {
    ctfN = ctf1;
  }

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
               const int & n, const Periodic & p,
               const BndGrid & co1, const BndGrid & coN) :
               nc_in(n), nc_tot(n+2*boil::BW),
               dummy_grid(false),
               period1(p), periodN(p), ctf1(co1), ctfN(coN) {

#if 1
   if(nc_in<boil::BW) {
     boil::aout<<"At least as many cells as the buffer width ("<<boil::BW
               <<") are required in each direction. Exiting."<<boil::endl;
     exit(0);
   }
#endif

/*----------------+
|  check cutoffs  |
+----------------*/
  if(p == Periodic::no()&&ctf1 == BndGrid::undefined()) {
    boil::aout<<"At least one cutoff must be specified for a non-periodic grid."
              <<" Exiting."<<boil::endl;
    exit(0);
  }
  if(p == Periodic::no()&&ctfN == BndGrid::undefined()) {
    ctfN = ctf1;
  }

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
               const Periodic & p, const BndGrid & co1, const BndGrid & coN) :
               nc_in(left.ncell()+right.ncell()), 
               nc_tot(nc_in+2*boil::BW),
               dummy_grid(false),
               period1(p), periodN(p), ctf1(co1), ctfN(coN) {
/*----------------+
|  check cutoffs  |
+----------------*/
  if(p == Periodic::no()&&ctf1 == BndGrid::undefined()) {
    boil::aout<<"At least one cutoff must be specified for a non-periodic grid."
              <<" Exiting."<<boil::endl;
    exit(0);
  }
  if(p == Periodic::no()&&ctfN == BndGrid::undefined()) {
    ctfN = ctf1;
  }

/*---------------------------+
|  creates non-uniform grid  |
+---------------------------*/

  allocate(); 

#if 0 /* this doesn't work with expanded buffers */
  /* just copty the first grid */
  for(int i=1; i<=left.nnode(); i++)
    x_node[i] = left.xn(i);

  /* append the second grid */ 
  for(int i=2; i<=right.nnode(); i++)
    x_node[i+left.ncell()] = x_node[left.nnode()] + right.xn(i) - right.xn(1);
#else
  /* copy first grid */
  for(int i=0; i<left.nnode(); i++)
    x_node[i+boil::BW] = left.xn(i+boil::BW);

  /* append the second grid */ 
  for(int i=1; i<right.nnode(); i++)
    x_node[i+left.nnode()+boil::BW-1] = x_node[left.nnode()+boil::BW-1]
                                      + right.xn(i+boil::BW)
                                      - right.xn(boil::BW);
#endif

  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D & left, const Grid1D & center, const Grid1D & right, 
               const Periodic & p, const BndGrid & co1, const BndGrid & coN) :
  nc_in(left.ncell()+center.ncell()+right.ncell()), 
  nc_tot(nc_in+2*boil::BW),
  dummy_grid(false),
  period1(p), periodN(p), ctf1(co1), ctfN(coN) {
/*----------------+
|  check cutoffs  |
+----------------*/
  if(p == Periodic::no()&&ctf1 == BndGrid::undefined()) {
    boil::aout<<"At least one cutoff must be specified for a non-periodic grid."
              <<" Exiting."<<boil::endl;
    exit(0);
  }
  if(p == Periodic::no()&&ctfN == BndGrid::undefined()) {
    ctfN = ctf1;
  }

/*---------------------------+
|  creates non-uniform grid  |
+---------------------------*/

  allocate();

#if 0 /* this doesn't work with expanded buffers */
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
#else
  /* copy first grid */
  for(int i=0; i<left.nnode(); i++)
    x_node[i+boil::BW] = left.xn(i+boil::BW);

  int first = boil::BW+left.nnode();
  int i0 = 1; /* one node is common */
  /* append the second grid */
  for(int i=i0; i<center.nnode(); i++)
    x_node[i+first-i0] = x_node[first-1]
                       + center.xn(i+boil::BW) - center.xn(boil::BW);

  first = boil::BW+left.nnode()+center.nnode()-i0; /* one node is common */
  /* append the third grid */
  for(int i=i0; i<right.nnode(); i++)
    x_node[i+first-i0] = x_node[first-1]
                       + right.xn(i+boil::BW) - right.xn(boil::BW);
#endif

  correct_boundaries();
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D & grid, 
               const Step   & step) // default step is Step(1)
  : period1(grid.periodic1()), 
    periodN(grid.periodicN()),
    dummy_grid(grid.is_dummy()),
    ctf1(grid.cutoff1()),
    ctfN(grid.cutoffN()) {
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

  if(grid.is_dummy()) {
    if(step.size() != 1) {
      boil::aout << "Critical! Can not coarsen grid! Stopping!" << boil::endl;
      exit(0);
    } else {
      nc_in = grid.ncell();
      nc_tot = nc_in + 2*boil::BW;

      real dx = grid.x_max() - grid.x_min();
      dummy_setup(dx);
    }
  } else {
    /* check if it is possible */
    if(grid.ncell() % step.size() !=0) {
      boil::aout << "Critical! Can not coarsen grid! Stopping!" << boil::endl;
      exit(0);
    }

    /* if yes, create (coarser) grid */
    nc_in = grid.ncell()/step.size();
    nc_tot = nc_in + 2*boil::BW;

#if 1
    if(nc_in<boil::BW) {
       boil::aout<<"At least as many cells as the buffer width ("<<boil::BW
                 <<") are required in each direction. Exiting."<<boil::endl;
       exit(0);
    }
#endif

    allocate(); 
  
    for(int i=0; i<=nc_in; i++)
      x_node[boil::BW + i] = grid.xn(boil::BW + i*step.size());

    correct_boundaries();
  }
}  

/******************************************************************************/
Grid1D::Grid1D(const Grid1D     & grid, 
               const Range<int> & cr,     // first and last cell inside
               const Step       & step) 
  : period1(grid.periodic1()), 
    periodN(grid.periodicN()),
    dummy_grid(grid.is_dummy()),
    ctf1(grid.cutoff1()),
    ctfN(grid.cutoffN()) {
/*------------------------------------------------------------------------+
|  copy constructor which gets a subgrid.                                 | 
|  the disadvantage is that it allocates and de-allocates memory :-(      |
+------------------------------------------------------------------------*/

  /* check if it is possible */
  if(grid.is_dummy()&&step.size()!=1) {
    boil::aout<<"Critical! No subgrids can be taken from a dummy grid! Stopping!"
              <<boil::endl;
    exit(0);
  }

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
  nc_tot = nc_in + 2*boil::BW;
  if(grid.is_dummy()) {
    real dx = grid.x_max() - grid.x_min();
    dummy_setup(dx);
  } else {
#if 1
    if(nc_in<boil::BW) {
       boil::aout<<"At least as many cells as the buffer width ("<<boil::BW
                 <<") are required in each direction. Exiting."<<boil::endl;
       exit(0);
    }
#endif

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
  }

  /* clean up */
  delete coarse;
}  

/******************************************************************************/
Grid1D::Grid1D(const real dx) :
  nc_in(1),
  nc_tot(1+2*boil::BW),
  dummy_grid(true),
  period1(Periodic::yes()),
  periodN(Periodic::yes()),
  ctf1(BndGrid::wall()),
  ctfN(BndGrid::wall()) {

  dummy_setup(dx);
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
  nc_tot = 0;
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
    if(i == boil::BW || i == nnode_b()-boil::BW) boil::aout << "- - - - - -" << boil::endl; 
    boil::aout << "i="<<i << "   xn="<<xn(i) << "   dxn="<<dxn(i) 
               << boil::endl;
  }

  if(periodicN()==true) boil::aout << ":::::::::::" << boil::endl;
  boil::aout << boil::endl;
  boil::aout << "-[ cell based dimensions ]-------------------------------" 
             << boil::endl;
  if(periodic1()==true) boil::aout << ":::::::::::" << boil::endl;
  for(int i=0; i<ncell_b(); i++) {
    if(i == boil::BW || i == ncell_b()-boil::BW) boil::aout << "- - - - - -" << boil::endl; 
    boil::aout << "i="<<i << "   xc="<<xc(i) << "   dxc="<<dxc(i) 
               << boil::endl;
  }
  if(periodicN()==true) boil::aout << ":::::::::::" << boil::endl;
}

