#ifndef DOMAIN_H
#define DOMAIN_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>

#include "../Grid/grid1d.h"
#include "../Body/body.h"
#include "../Body/Empty/empty.h"
#include "../Ravioli/periodic.h"
#include "../Ravioli/bndgrid.h"
#include "../Ravioli/dir.h"
#include "../Ravioli/comp.h"
#include "../Ravioli/decompose.h"
#include "../Global/global_swap.h"

//////////////
//          //
//  Domain  //
//          //
//////////////
class Domain {
  public:
    Domain(const Grid1D & ogx, const Grid1D & ogy, const Grid1D & ogz,
           const std::string n="domain", const Decompose dec=Decompose::xyz());

    Domain(const Grid1D & ogx, const Grid1D & ogy, const Grid1D & ogz,
           Body * b, /* it will change, that is why it is pointer */ 
           const std::string n="domain", const Decompose dec=Decompose::xyz());

    ~Domain(){}

    int  level() const {return lev;}

    /* local number of cells */
    int  ni()    const {return grid_x_local->ncell_b();}
    int  nj()    const {return grid_y_local->ncell_b();}
    int  nk()    const {return grid_z_local->ncell_b();}

    /* global number of cells */
    int  gi() const {return grid_x_original->ncell_b();}
    int  gj() const {return grid_y_original->ncell_b();}
    int  gk() const {return grid_z_original->ncell_b();}


    /*  global extents of the domain */
    const Range<int> & cxg() const {return cr_x;}
    const Range<int> & cyg() const {return cr_y;}
    const Range<int> & czg() const {return cr_z;}

    /* cell dimensions */
    real dxc(const int i) const {return grid_x_local->dxc(i);}
    real dxw(const int i) const {return grid_x_local->dxn(i);}
    real dxe(const int i) const {return grid_x_local->dxn(i+1);}
    real dyc(const int j) const {return grid_y_local->dxc(j);}
    real dys(const int j) const {return grid_y_local->dxn(j);}
    real dyn(const int j) const {return grid_y_local->dxn(j+1);}
    real dzc(const int k) const {return grid_z_local->dxc(k);}
    real dzb(const int k) const {return grid_z_local->dxn(k);}
    real dzt(const int k) const {return grid_z_local->dxn(k+1);}

    /* cell centre coordinates */
    real xc(const int i) const {return grid_x_local->xc(i);}
    real yc(const int j) const {return grid_y_local->xc(j);}
    real zc(const int k) const {return grid_z_local->xc(k);}

    /* node coordinates */
    real xn(const int i) const {return grid_x_local->xn(i);}
    real yn(const int j) const {return grid_y_local->xn(j);}
    real zn(const int k) const {return grid_z_local->xn(k);}

    /* cell centre coordinates */
    real xc_global(const int I) const {return grid_x_original->xc(I);}
    real yc_global(const int J) const {return grid_y_original->xc(J);}
    real zc_global(const int K) const {return grid_z_original->xc(K);}

    /* node coordinates */
    real xn_global(const int I) const {return grid_x_original->xn(I);}
    real yn_global(const int J) const {return grid_y_original->xn(J);}
    real zn_global(const int K) const {return grid_z_original->xn(K);}

    /* global coordinates extents */
    real global_min_x() const {return grid_x_original->x_min();}
    real global_min_y() const {return grid_y_original->x_min();}
    real global_min_z() const {return grid_z_original->x_min();}
    real global_max_x() const {return grid_x_original->x_max();}
    real global_max_y() const {return grid_y_original->x_max();}
    real global_max_z() const {return grid_z_original->x_max();}

    /* carefull: these return global logical coordinates */
    int I(const real x) const;
    int J(const real y) const;
    int K(const real z) const;

    /* cell surfaces */
    real dSx(const int i, const int j, const int k) const 
     {return dyc(j) * dzc(k);}
    real dSy(const int i, const int j, const int k) const
     {return dxc(i) * dzc(k);}
    real dSz(const int i, const int j, const int k) const
     {return dxc(i) * dyc(j);}

    /* cell volume */
    real dV(const int i, const int j, const int k) const 
     {return dxc(i) * dyc(j) * dzc(k);}

    bool  period(const int i) const {return per[i];}
    bool  cutoff(const int i, const int j) const {return ctf[i][j];}
    bool  bnd_symmetry(const Dir d) const;
    const Domain * coarser() const {return crsr;}
    
    const Body & ibody() const {return * body;}

    /* these functions check the global cell range (excluding buffers) */
    bool contains_I(int I) const {return cr_x.contains(I);}
    bool contains_J(int J) const {return cr_y.contains(J);}
    bool contains_K(int K) const {return cr_z.contains(K);}
    bool contains_IJK(int I, int J, int K) const {
      return contains_I(I) && contains_J(J) && contains_K(K);
    } 
    /* these functions check the local cell range (excluding buffers) */
    bool contains_i(int i) const {return i>0 && i<ni()-1;}
    bool contains_j(int j) const {return j>0 && j<nj()-1;}
    bool contains_k(int k) const {return k>0 && k<nk()-1;}
    bool contains_ijk(int i, int j, int k) const {
      return contains_i(i) && contains_j(j) && contains_k(k);
    } 
    /* these functions check the local coordinate range (excluding buffers) */
    bool contains_x(real x) const { return x>=xn(1) && x<=xn(ni()-1);}
    bool contains_y(real y) const { return y>=yn(1) && y<=yn(nj()-1);}
    bool contains_z(real z) const { return z>=zn(1) && z<=zn(nk()-1);}
    bool contains_xyz(real x, real y, real z) const {
      return contains_x(x) && contains_y(y) && contains_z(z);
    } 

    void locals(int * i, int * j, int * k) const; 
    int local_i(int I) const; 
    int local_j(int J) const; 
    int local_k(int K) const; 
    void globals(int * i, int * j, int * k) const; 
    int global_I(int i) const; 
    int global_J(int j) const; 
    int global_K(int k) const; 

    /* are these needed? */
    int neighbour(const int n) const {return neighbours[n];}
    int coord(const Comp & i)  const {return coords[~i];}
    int dim(const Comp & i)    const {return dims[~i];}

    /* computes and prints grid statistics */
    void statistics(Body * b = NULL);

    real dxyz_min() const {return min_dxyz;}
    real dxyz_max() const {return max_dxyz;}
    real dV_min()   const {return min_dV;}
    real dV_max()   const {return max_dV;}

  private:
    Domain(const Grid1D * ogx, const Grid1D * ogy, const Grid1D * ogz,
           const Grid1D * lgx, const Grid1D * lgy, const Grid1D * lgz,
           int * dms, int * crds, int * nghbrs,
           const int l, 
           const Range<int> & crx, const Range<int> & cry, const Range<int> & crz,
           const std::string & n);
    void distribute(const int nproc, const int dim, int * dis, int * res);
    void factor(int n, int * factor, int * number) const;
    void init(const Decompose & dec);
    void decompose(const int i, const int g, Range<int> * cr) const;

    const Grid1D * grid_x_original, * grid_y_original, * grid_z_original;
    const Grid1D * grid_x_local,    * grid_y_local,    * grid_z_local;    

    Body * body; /* immersed body */

    const int lev;               /* refinement level */
    Range<int> cr_x, cr_y, cr_z; /* cell range - for domain decomposition */
    bool  per[3];
    bool  ctf[3][2];
    const Domain * crsr;
    const Domain * coarsen() const;
    void  setup(const Decompose & dec);

    /* for domain decomposition */
    int * neighbours; /* neighbouring subdomains */
    int * coords;
    int * dims;       /* number of sub-domains in each direction */

    /* name */
    std::string name;

    /* for grid statistics */
    real min_dx, max_dx, min_dy, max_dy, min_dz, max_dz;
    real min_dxyz, max_dxyz;
    real min_dV, max_dV;
    real max_ar; /* aspect ratio */
    bool imins,imaxs,jmins,jmaxs,kmins,kmaxs;
};	

#endif
