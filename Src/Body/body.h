#ifndef BODY_H
#define BODY_H

#include "../Parallel/mpi_macros.h"
#include <fstream>
#include <vector>
#include <cfloat>
#include <algorithm>

#include "../Parallel/Out/out.h"
#include "../Parallel/Out/print.h"
#include "../Global/global_precision.h"
#include "../Global/global_approx.h"
#include "../Global/global_minmax.h"
#include "../Ravioli/range.h"
#include "../Ravioli/comp.h"
#include "../Field/Vector/vector_browsing.h" 
#include "../Timer/timer.h"

#include "anglenode.h"
#include "cutcell.h"
#include "polygon.h"
#include "plane.h"

#include "body_nonmember.h"

class Domain;
class Scalar;
class Vector;
class VectorBool;

////////////
//        //
//  Body  //
//        //
////////////
class Body {
  public:
    Body(const std::string name);
    Body() {sca=NULL; vec=NULL;}
    //! Number of Polygons (cells).
    int npolys() const {return polys.size();}

    //! Total (over all domains) number of Polygons (cells).
    int tpolys() const {int t=npolys(); boil::cart.sum_int(&t); return t;}

    //! Total number of nodes (summ of all Polygons' (cells') nodes).
    int nnodes() const {
      int t=0; 
      for(int c=0; c<npolys(); c++) t+=nnodes(c);
      return t;}                

    //! Number of nodes of polygon "c". 
    int nnodes(const int c) const {return polys[c].nnodes();}

    //! Polygon's "c" normal component "i". 
    real n(const int m, const int c) const {return polys[c].n(m);}

    //! Polygon's "c" surface components. 
    real dSx(const int c) const {return polys[c].area_x();}
    real dSy(const int c) const {return polys[c].area_y();}
    real dSz(const int c) const {return polys[c].area_z();}
    real dSx(const int i, const int j, const int k) const;
    real dSy(const int i, const int j, const int k) const;
    real dSz(const int i, const int j, const int k) const;

    //! Polygon's "c" surface total area. 
    real dS(const int c) const {
      return sqrt( dSx(c)*dSx(c) + dSy(c)*dSy(c) + dSz(c)*dSz(c) );
    }
    real dS(const int i, const int j, const int k) const {
      return sqrt(   dSx(i,j,k)*dSx(i,j,k) 
                   + dSy(i,j,k)*dSy(i,j,k) 
                   + dSz(i,j,k)*dSz(i,j,k) );
    }

    const Polygon & operator [] (const int c) const {return polys[c];}

    void cut(const Domain & dom); 

    void cut(const Scalar & phi); 

    //! For communication with Equation
    int nccells() const;
    int ncall() const;
    int nccells(const Comp & m) const;
    void ijk(int cc, int * i, int * j, int * k) const;  
    void ijk(Comp & m, int cc, int * i, int * j, int * k) const;  

    virtual bool cut(int i, int j, int k) const;
    bool cut(int cc) const;
    virtual bool off(int i, int j, int k) const;
    virtual bool on (int i, int j, int k) const;

    /* only for pressure equation */
    virtual bool cut_p(int i, int j, int k) const;
    virtual bool off_p(int i, int j, int k) const;
    virtual bool on_p (int i, int j, int k) const;

    virtual bool cut(const Comp & m, int i, int j, int k) const;
    bool cut(const Comp & m, int cc) const;
    virtual bool off(const Comp & m, int i, int j, int k) const;
    virtual bool on (const Comp & m, int i, int j, int k) const;

    real fSw(const int cc) const; 
    real fSe(const int cc) const; 
    real fSs(const int cc) const; 
    real fSn(const int cc) const; 
    real fSb(const int cc) const; 
    real fSt(const int cc) const; 

    real fSw(int i, int j, int k) const; 
    real fSe(int i, int j, int k) const; 
    real fSs(int i, int j, int k) const; 
    real fSn(int i, int j, int k) const; 
    real fSb(int i, int j, int k) const; 
    real fSt(int i, int j, int k) const; 

    real fSw(const Comp & m, const int cc) const; 
    real fSe(const Comp & m, const int cc) const; 
    real fSs(const Comp & m, const int cc) const; 
    real fSn(const Comp & m, const int cc) const; 
    real fSb(const Comp & m, const int cc) const; 
    real fSt(const Comp & m, const int cc) const; 

    virtual real fSw(const Comp & m, int i, int j, int k) const; 
    virtual real fSe(const Comp & m, int i, int j, int k) const; 
    virtual real fSs(const Comp & m, int i, int j, int k) const; 
    virtual real fSn(const Comp & m, int i, int j, int k) const; 
    virtual real fSb(const Comp & m, int i, int j, int k) const; 
    virtual real fSt(const Comp & m, int i, int j, int k) const; 

    real fdxw(const int cc) const; 
    real fdxe(const int cc) const; 
    real fdys(const int cc) const; 
    real fdyn(const int cc) const; 
    real fdzb(const int cc) const; 
    real fdzt(const int cc) const; 

    real fdxw(int i, int j, int k) const;
    real fdxe(int i, int j, int k) const;
    real fdys(int i, int j, int k) const;
    real fdyn(int i, int j, int k) const;
    real fdzb(int i, int j, int k) const;
    real fdzt(int i, int j, int k) const;

    real fdxw(const Comp & m, const int cc) const; 
    real fdxe(const Comp & m, const int cc) const; 
    real fdys(const Comp & m, const int cc) const; 
    real fdyn(const Comp & m, const int cc) const; 
    real fdzb(const Comp & m, const int cc) const; 
    real fdzt(const Comp & m, const int cc) const; 

    real fPmmm(const int cc) const;
    real fPpmm(const int cc) const;
    real fPmpm(const int cc) const;
    real fPppm(const int cc) const;
    real fPmmp(const int cc) const;
    real fPpmp(const int cc) const;
    real fPmpp(const int cc) const;
    real fPppp(const int cc) const;

    real fPmmm(int i, int j, int k) const;
    real fPpmm(int i, int j, int k) const;
    real fPmpm(int i, int j, int k) const;
    real fPppm(int i, int j, int k) const;
    real fPmmp(int i, int j, int k) const;
    real fPpmp(int i, int j, int k) const;
    real fPmpp(int i, int j, int k) const;
    real fPppp(int i, int j, int k) const;

    real fE000(const int cc) const;
    real fE010(const int cc) const;
    real fE001(const int cc) const;
    real fE011(const int cc) const;
    real fE100(const int cc) const;
    real fE110(const int cc) const;
    real fE101(const int cc) const;
    real fE111(const int cc) const;
    real fE200(const int cc) const;
    real fE210(const int cc) const;
    real fE201(const int cc) const;
    real fE211(const int cc) const;

    real fE000(int i, int j, int k) const;
    real fE010(int i, int j, int k) const;
    real fE001(int i, int j, int k) const;
    real fE011(int i, int j, int k) const;
    real fE100(int i, int j, int k) const;
    real fE110(int i, int j, int k) const;
    real fE101(int i, int j, int k) const;
    real fE111(int i, int j, int k) const;
    real fE200(int i, int j, int k) const;
    real fE210(int i, int j, int k) const;
    real fE201(int i, int j, int k) const;
    real fE211(int i, int j, int k) const;

    real fV(int i, int j, int k) const; 
    real fV(const int cc) const; 

    real fV(const Comp & m, int i, int j, int k) const;
    real fV(const Comp & m, const int cc) const;

    real dist(const int i, const int j, const int k) const;
    real nwx(const int i, const int j, const int k) const;
    real nwy(const int i, const int j, const int k) const;
    real nwz(const int i, const int j, const int k) const;

    int iacpt(const int cc) const;
    int jacpt(const int cc) const;
    int kacpt(const int cc) const;
    int iacpt(const int i, const int j, const int k) const;
    int jacpt(const int i, const int j, const int k) const;
    int kacpt(const int i, const int j, const int k) const;

    Scalar * bdist;
    Scalar * ux, * uy, * uz;

  private:

    bool cross_seg_x(int index,
                     const Range<real> x_seg, const real y, const real z,
                     real * f,
                     int  * cp) const;
    bool cross_seg_y(int index,
                     const real x, const Range<real> y_seg, const real z,
                     real * f,
                     int  * cp) const;
    bool cross_seg_z(int index,
                     const real x, const real y, const Range<real> z_seg,
                     real * f,
                     int  * cp) const;

    CutCell * cut_cell(int index,
                       int i, int j, int k,
                       const real x0, const real y0, const real z0,
                       const real dx, const real dy, const real dz,
                       const real xc, const real yc, const real zc,
                       const real xw, const real xe, 
                       const real ys, const real yn,
                       const real zb, const real zt,
                       int m,
                       Body * new_faces = NULL);

    Scalar * sca; /* helping array */
    Vector * vec; /* helping vector */
    VectorBool * vecoff;
    Scalar * stmp, * dflag; /* helping array */
    std::vector<int ***> index;

    std::vector<Polygon>                polys; 
    std::vector< std::vector<int> >     polytags; /* poly's for each cell */
    std::vector< std::vector<CutCell> > cells;

    int nccells_in[4], nccells_bd[4];
    int ncall_;

    /* core of cutting */
    void cut_init(const Domain & dom);
    void cut_center(const Domain & dom);
    void cut_stgd(const Domain & dom);
    void set_vecoff();

    /* tool for cutting */
    void cut_degen(real x[], real y[], real z[], int & i,
                   int ic, int jc, int kc);

    /* set volume (fV) */
    void vol_center(const Domain & dom, const Body & nf);
    void vol_stgd(const Domain & dom, const Body & b1, const Body & b2,
                   const Body & b3);

    /* set area (fS) */
    void fluid_area(CutCell *c, int i[][2][2], bool b[][2][2],
                    real x1[][2], real y1[][2], real z1[][2],
                    real x2[], real y2[], real z2[],
                    real dx,   real xy,   real dz );

    /* set edge ratio (fE) */
    void cut_edge_ratio(CutCell * ccell,Plane * p, 
           const int n_in_fluid[][2][2],
           real x[], real y[], real z[]);

    /* check face-match */
    void match_face(int index, int i, int j, int k
                      ,real x[], real y[], real z[], int face_match[][2]
                      ,real face_match_norm[][2][3], int & nf);

    /* functions and variables for floodfill */
    void body_fills(real new_color, int i, int j, int k);
    void body_fillv(real new_color, int i, int j, int k, const Comp m);
    bool pop(int * i, int * j, int * k);
    void push(int i, int j, int k);

    /* distance function */
    void distfunc(const Domain & d);
    void nwall(const Domain & d);
    void insert_bc_dist(const Scalar * sca);

    /* cell merging */
    void acceptor();

    /* boundary cell treatment */
    void bd_center(const Domain & dom, Scalar * sca);

    /* check duplication */
    int cut_duplic(Plane * p, const real x[], const real y[],
                   const real z[], const int n_in_fluid[][2][2]) const;

    void normalize(real & r1, real & r2, real & r3);

    int stack_size;
    int ** stack;
    int stack_pointer;
    int sifl[4],eifl[4],sjfl[4],ejfl[4],skfl[4],ekfl[4];
    real tol, dxmin;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: body.h,v 1.45 2015/08/17 10:33:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
