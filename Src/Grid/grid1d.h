#ifndef GRID1D_H
#define GRID1D_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include "../Global/global_checking.h"
#include "../Global/global_endl.h"
#include "../Global/global_precision.h"
#include "../Ravioli/periodic.h"
#include "../Ravioli/range.h"
#include "../Parallel/Out/out.h"
#include "../Parallel/Out/print.h"
#include "../Board/board.h"

#include "step.h"

//////////////
//          //
//  Grid1D  //
//          //
//////////////
class Grid1D {
  public:
    Grid1D(const Range<real> &  xr, 
           const Range<real> & dxr,
           const int  & gx,
           const Periodic & p);

    Grid1D(const Range<real> &  xr, 
           const int  & gx,
           const Periodic & p);

    Grid1D(const Grid1D & left,
           const Grid1D & right,
           const Periodic & p);

    Grid1D(const Grid1D & left,
           const Grid1D & center, 
           const Grid1D & right,
           const Periodic & p);

    /* constructor which creates a coarser grid */
    Grid1D(const Grid1D & grid, const Step & step = Step(1));

    /* constructor which creates a subgrid */
    Grid1D(const Grid1D     & grid, 
           const Range<int> & cr, 
           const Step       & step = Step(1));
    ~Grid1D();

    /* are these used? */
    const Periodic periodic1() const {return period1;}
    const Periodic periodicN() const {return periodN;}

    int ncell()   const {return N;}
    int nnode()   const {return N+1;}
    int ncell_b() const {return N+2;}
    int nnode_b() const {return N+3;}

    real xn (const int i) const {return  x_node[i];}
    real xc (const int i) const {return  x_cell[i];}
    real dxn(const int i) const {return dx_node[i];}
    real dxc(const int i) const {return dx_cell[i];}
    
    void print() const;
    void plot(const char *) const;

    real x_min() const {return x_node[1];}
    real x_max() const {return x_node[N+1];}

    bool contains(real x) const {
      if( x >= x_min() && x <= x_max() ) return true;
      return false;
    }

  private:
    int    N;       /* number of cells inside!!! */
    real *  x_node;
    real *  x_cell;
    real * dx_node;
    real * dx_cell;
    Periodic period1, periodN;
    const    Grid1D * coarsen() const;
    void     allocate();
    void     distribute_nodes_inside(const real & x1, const real & xn,
                                     const real & D1, const real & DN,
                                     const int  & gx);
    void     correct_boundaries();
};

#endif
