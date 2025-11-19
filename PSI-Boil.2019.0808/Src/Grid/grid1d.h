#ifndef GRID1D_H
#define GRID1D_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include "../Global/global_checking.h"
#include "../Global/global_endl.h"
#include "../Global/global_precision.h"
#include "../Ravioli/bndgrid.h"
#include "../Ravioli/periodic.h"
#include "../Ravioli/range.h"
#include "../Parallel/Out/out.h"
#include "../Parallel/Out/print.h"

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
           const Periodic & p,
           const BndGrid & co1 = BndGrid::wall(),
           const BndGrid & coN = BndGrid::undefined());

    Grid1D(const Range<real> &  xr, 
           const int  & gx,
           const Periodic & p,
           const BndGrid & co1 = BndGrid::wall(),
           const BndGrid & coN = BndGrid::undefined());

    Grid1D(const Grid1D & left,
           const Grid1D & right,
           const Periodic & p,
           const BndGrid & co1 = BndGrid::wall(),
           const BndGrid & coN = BndGrid::undefined());

    Grid1D(const Grid1D & left,
           const Grid1D & center, 
           const Grid1D & right,
           const Periodic & p,
           const BndGrid & co1 = BndGrid::wall(),
           const BndGrid & coN = BndGrid::undefined());

    /* constructor to create a dummy grid (for lower dim simulations) */
    explicit Grid1D(const real dx);

    /* constructor which creates a coarser grid */
    Grid1D(const Grid1D & grid, const Step & step = Step(1));

    /* constructor which creates a subgrid */
    Grid1D(const Grid1D     & grid, 
           const Range<int> & cr, 
           const Step       & step = Step(1));
    ~Grid1D();

    const Periodic periodic1() const {return period1;}
    const Periodic periodicN() const {return periodN;}

    const BndGrid bndgrid1() const {return ctf1;}
    const BndGrid bndgridN() const {return ctfN;}

    /* number of cells and nodes */
    int ncell()   const {return nc_in;}
    int nnode()   const {return nc_in + 1;}

    /* number of cells and nodes including boundaries (or buffers) */
    int ncell_b() const {return nc_tot;}
    int nnode_b() const {return nc_tot + 1;}

    real xn (const int i) const {return  x_node[i];}
    real xc (const int i) const {return  x_cell[i];}
    real dxn(const int i) const {return dx_node[i];}
    real dxc(const int i) const {return dx_cell[i];}

    void print() const;

    real x_min() const {
      return x_node[boil::BW];
    }
    real x_max() const {
      return x_node[nc_in + boil::BW];
    }

    real lx() const {
      return x_max() - x_min();
    }

    bool contains(real x) const {
      if( x >= x_min() && x <= x_max() ) return true;
      return false;
    }

    bool is_dummy() const { return dummy_grid; }

  private:
    int     nc_in;       /* number of cells inside!!! */
    int     nc_tot;      /* number of cells in total  */
    bool dummy_grid;     /* is this a dummy grid? */
    real *  x_node;
    real *  x_cell;
    real * dx_node;
    real * dx_cell;
    Periodic period1, periodN;
    BndGrid ctf1, ctfN;
    const    Grid1D * coarsen() const;
    void     allocate();
    void     distribute_nodes_inside(const real & x1, const real & xn,
                                     const real & D1, const real & DN);
    void     dummy_setup(const real dx);
    void     correct_boundaries();
};

#endif
