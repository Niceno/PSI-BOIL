#ifndef VECTORBOOL_H
#define VECTORBOOL_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include "../../Ravioli/comp.h"
#include "../../Ravioli/sign.h"
#include "../../Parallel/communicator.h"
#include "../../Domain/domain.h"
#include "../../Boundary/bndcnd.h"
#include "../Vector/vector_browsing.h"
#include "../../Global/global_malloc.h"
#include "../../Global/global_minmax.h"
#include "../../Global/global_approx.h"
#include "../ScalarBool/scalarbool.h"

//#include "../../Timer/timer.h"

////////////////////////////////
//                            //
//  ScalarP (Scalar pointer)  //
//                            //
////////////////////////////////
class ScalarBoolP {
  public:
    ScalarBoolP() {sp = new ScalarBool[3];}
    ScalarBool * sp;
    ScalarBool & operator [] (const Comp & m) const {return sp[~m];}
};

//////////////
//          //
//  Vector  //
//          //
//////////////
class VectorBool {
  public:
    /* global constructor */
    VectorBool(const Domain & d);
    VectorBool(const VectorBool * v); // creates an alias

    /* local constructors and functions */
    VectorBool() : aliav(true) {};
    void allocate(const int ni, const int nj, const int nk);
    ~VectorBool();
    
    int ni() const {return dom->ni();}
    int nj() const {return dom->nj();}
    int nk() const {return dom->nk();}

    int ni(const Comp & m) const {return vec[m].n_x;}
    int nj(const Comp & m) const {return vec[m].n_y;}
    int nk(const Comp & m) const {return vec[m].n_z;}
    int ox(const Comp & m) const {return vec[m].o_x;}
    int oy(const Comp & m) const {return vec[m].o_y;}
    int oz(const Comp & m) const {return vec[m].o_z;}
    int si(const Comp & m) const {return vec[m].s_x;}
    int sj(const Comp & m) const {return vec[m].s_y;}
    int sk(const Comp & m) const {return vec[m].s_z;}
    int ei(const Comp & m) const {return vec[m].e_x;}
    int ej(const Comp & m) const {return vec[m].e_y;}
    int ek(const Comp & m) const {return vec[m].e_z;}

    int & si(const Comp & m) {return vec[m].s_x;}
    int & sj(const Comp & m) {return vec[m].s_y;}
    int & sk(const Comp & m) {return vec[m].s_z;}
    int & ei(const Comp & m) {return vec[m].e_x;}
    int & ej(const Comp & m) {return vec[m].e_y;}
    int & ek(const Comp & m) {return vec[m].e_z;}

    /* scalar node coordinates (needed for plotting) */
    real xn(const int i) const {return dom->xn(i);}
    real yn(const int j) const {return dom->yn(j);}
    real zn(const int k) const {return dom->zn(k);}

    /* scalar cell dimensions (needed for velocity gradients) */
    real dxc(const int i) const {return dom->dxc(i);}
    real dyc(const int j) const {return dom->dyc(j);}
    real dzc(const int k) const {return dom->dzc(k);}

    real dxw(const int i) const {return dom->dxw(i);}
    real dxe(const int i) const {return dom->dxe(i);}
    real dys(const int j) const {return dom->dys(j);}
    real dyn(const int j) const {return dom->dyn(j);}
    real dzb(const int k) const {return dom->dzb(k);}
    real dzt(const int k) const {return dom->dzt(k);}

    /* vector node coordinates */
    real xn(const Comp & m, const int i) const {return (this->*pnt_xn[~m])(i);}
    real yn(const Comp & m, const int j) const {return (this->*pnt_yn[~m])(j);}
    real zn(const Comp & m, const int k) const {return (this->*pnt_zn[~m])(k);}

    /* vector cell dimensions */
    real dxc(const Comp & m, const int i) const {return (this->*pnt_dxc[~m])(i);}
    real dyc(const Comp & m, const int j) const {return (this->*pnt_dyc[~m])(j);}
    real dzc(const Comp & m, const int k) const {return (this->*pnt_dzc[~m])(k);}

    real dxw(const Comp & m, const int i) const {return (this->*pnt_dxw[~m])(i);}
    real dxe(const Comp & m, const int i) const {return (this->*pnt_dxe[~m])(i);}
    real dys(const Comp & m, const int j) const {return (this->*pnt_dys[~m])(j);}
    real dyn(const Comp & m, const int j) const {return (this->*pnt_dyn[~m])(j);}
    real dzb(const Comp & m, const int k) const {return (this->*pnt_dzb[~m])(k);}
    real dzt(const Comp & m, const int k) const {return (this->*pnt_dzt[~m])(k);}

    /* vector cell coordinates */
    real  xc(const Comp & m, const int i) const {return (this->*pnt_xc[~m])(i);}
    real  yc(const Comp & m, const int j) const {return (this->*pnt_yc[~m])(j);}
    real  zc(const Comp & m, const int k) const {return (this->*pnt_zc[~m])(k);}

    /* staggered cell-face areas */
#if 0
    real dSx(const Comp & m, const int i, const int j, const int k) const
     {return dyc(m,j) * dzc(m,k);}
    real dSy(const Comp & m, const int i, const int j, const int k) const
     {return dxc(m,i) * dzc(m,k);}
    real dSz(const Comp & m, const int i, const int j, const int k) const
     {return dxc(m,i) * dyc(m,j);}
#else
    real dSx(const Comp & m, const int i, const int j, const int k) const
     {return (this->*pnt_dSx[~m])(i,j,k);}
    real dSy(const Comp & m, const int i, const int j, const int k) const
     {return (this->*pnt_dSy[~m])(i,j,k);}
    real dSz(const Comp & m, const int i, const int j, const int k) const
     {return (this->*pnt_dSz[~m])(i,j,k);}
#endif

    real dSx(const Comp & m, const Sign sig, const int i, const int j, const int k) const
      {return (this->*pnt_dSx_dir[~m])(sig,i,j,k);}
    real dSy(const Comp & m, const Sign sig, const int i, const int j, const int k) const
      {return (this->*pnt_dSy_dir[~m])(sig,i,j,k);}
    real dSz(const Comp & m, const Sign sig, const int i, const int j, const int k) const
      {return (this->*pnt_dSz_dir[~m])(sig,i,j,k);}

    /* staggered cell volume */
#if 0
    real dV(const Comp & m, const int i, const int j, const int k) const
     {return dxc(m,i) * dyc(m,j) * dzc(m,k);}
#else
    real dV(const Comp & m, const int i, const int j, const int k) const
     {return (this->*pnt_dV[~m])(i,j,k);}
#endif

    int i(const Comp & m, const real x) const;
    int j(const Comp & m, const real y) const;
    int k(const Comp & m, const real z) const;
    int im(const Comp & m, const real x, const real t) const;
    int jm(const Comp & m, const real y, const real t) const;
    int km(const Comp & m, const real z, const real t) const;
    int ip(const Comp & m, const real x, const real t) const;
    int jp(const Comp & m, const real y, const real t) const;
    int kp(const Comp & m, const real z, const real t) const;
    int aim(const Comp & m, const real x, const real t) const;
    int ajm(const Comp & m, const real y, const real t) const;
    int akm(const Comp & m, const real z, const real t) const;
    int aip(const Comp & m, const real x, const real t) const;
    int ajp(const Comp & m, const real y, const real t) const;
    int akp(const Comp & m, const real z, const real t) const;

    void exchange    (const Comp comp = Comp::undefined(),
                      const int  dir  = -1);
    void exchange    (const int * i,
                      const Comp comp = Comp::undefined(),
                      const int  dir  = -1);
    void exchange_all(const Comp comp = Comp::undefined(),
                      const int  dir  = -1);

    bool *** operator [] (const Comp & m) const {return vec[m].val;}
    ScalarBool & operator () (const Comp & m) const {return vec[m];}


    /* mathematical operators */
    /* 
       = */
    const VectorBool & operator = (const bool * d)
     {for_m(m) 
        for_mijk(m,i,j,k) vec[m][i][j][k]=d[~m]; return *this;}
  
    void save(const char *, const int);
    void load(const char *, const int);
    void rm  (const char *, const int);
    void save(std::ofstream &);
    void load(std::ifstream &);

    BndCnd & bc(const Comp & m) const {return * (vec[m].bndcnd);}
    void bc_add(const BndCnd & bc) {for_m(m) vec[m].bndcnd->add(bc);}
    const Domain * domain() const {return dom;}

  private:
    void     coordinate();

    /* function pointers */
    real (VectorBool:: *pnt_xc[3])(const int i) const;
    real (VectorBool:: *pnt_yc[3])(const int j) const;
    real (VectorBool:: *pnt_zc[3])(const int k) const;

    real (VectorBool:: *pnt_xn[3])(const int i) const;
    real (VectorBool:: *pnt_yn[3])(const int j) const;
    real (VectorBool:: *pnt_zn[3])(const int k) const;

    real (VectorBool::*pnt_dxc[3])(const int i) const;
    real (VectorBool::*pnt_dyc[3])(const int j) const;
    real (VectorBool::*pnt_dzc[3])(const int k) const;

    real (VectorBool::*pnt_dxw[3])(const int i) const;
    real (VectorBool::*pnt_dxe[3])(const int i) const;

    real (VectorBool::*pnt_dys[3])(const int j) const;
    real (VectorBool::*pnt_dyn[3])(const int j) const;

    real (VectorBool::*pnt_dzb[3])(const int k) const;
    real (VectorBool::*pnt_dzt[3])(const int k) const;

    real (VectorBool::*pnt_dSx[3])(const int i, const int j, const int k) const;
    real (VectorBool::*pnt_dSy[3])(const int i, const int j, const int k) const;
    real (VectorBool::*pnt_dSz[3])(const int i, const int j, const int k) const;

    real (VectorBool::*pnt_dSx_dir[3])(const Sign sig, const int i, const int j, const int k) const;
    real (VectorBool::*pnt_dSy_dir[3])(const Sign sig, const int i, const int j, const int k) const;
    real (VectorBool::*pnt_dSz_dir[3])(const Sign sig, const int i, const int j, const int k) const;

    real (VectorBool::*pnt_dV[3])(const int i, const int j, const int k) const;

    /* these will be pointed to by above pointers */
    real  xc_nrm      (const int i) const {return dom-> xc(i);}
    real  xc_staggered(const int i) const {return dom-> xn(i);}
    real  yc_nrm      (const int j) const {return dom-> yc(j);}
    real  yc_staggered(const int j) const {return dom-> yn(j);}
    real  zc_nrm      (const int k) const {return dom-> zc(k);}
    real  zc_staggered(const int k) const {return dom-> zn(k);}

    real  xn_nrm      (const int i) const {return dom-> xn(i);}
    real  xn_staggered(const int i) const {return dom-> xc(i-1);}
    real  yn_nrm      (const int j) const {return dom-> yn(j);}
    real  yn_staggered(const int j) const {return dom-> yc(j-1);}
    real  zn_nrm      (const int k) const {return dom-> zn(k);}
    real  zn_staggered(const int k) const {return dom-> zc(k-1);}

    real dxc_nrm      (const int i) const {return dom->dxc(i);}
    real dxc_staggered(const int i) const {return dom->dxw(i);}
    real dyc_nrm      (const int j) const {return dom->dyc(j);}
    real dyc_staggered(const int j) const {return dom->dys(j);}
    real dzc_nrm      (const int k) const {return dom->dzc(k);}
    real dzc_staggered(const int k) const {return dom->dzb(k);}

    real dxw_nrm      (const int i) const {return dom->dxw(i);}
    real dxw_staggered(const int i) const {return dom->dxc(i-1);}
    real dxe_nrm      (const int i) const {return dom->dxe(i);}
    real dxe_staggered(const int i) const {return dom->dxc(i);}

    real dys_nrm      (const int j) const {return dom->dys(j);}
    real dys_staggered(const int j) const {return dom->dyc(j-1);}
    real dyn_nrm      (const int j) const {return dom->dyn(j);}
    real dyn_staggered(const int j) const {return dom->dyc(j);}

    real dzb_nrm      (const int k) const {return dom->dzb(k);}
    real dzb_staggered(const int k) const {return dom->dzc(k-1);}
    real dzt_nrm      (const int k) const {return dom->dzt(k);}
    real dzt_staggered(const int k) const {return dom->dzc(k);}

    real dSx_xstaggered(const int i, const int j, const int k) const {return dom->dSx_xstag(i,j,k);}
    real dSx_ystaggered(const int i, const int j, const int k) const {return dom->dSx_ystag(i,j,k);}
    real dSx_zstaggered(const int i, const int j, const int k) const {return dom->dSx_zstag(i,j,k);}

    real dSy_xstaggered(const int i, const int j, const int k) const {return dom->dSy_xstag(i,j,k);}
    real dSy_ystaggered(const int i, const int j, const int k) const {return dom->dSy_ystag(i,j,k);}
    real dSy_zstaggered(const int i, const int j, const int k) const {return dom->dSy_zstag(i,j,k);}

    real dSz_xstaggered(const int i, const int j, const int k) const {return dom->dSz_xstag(i,j,k);}
    real dSz_ystaggered(const int i, const int j, const int k) const {return dom->dSz_ystag(i,j,k);}
    real dSz_zstaggered(const int i, const int j, const int k) const {return dom->dSz_zstag(i,j,k);}

    real dSx_dir_xstaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSx_xstag(sig,i,j,k);}
    real dSx_dir_ystaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSx_ystag(sig,i,j,k);}
    real dSx_dir_zstaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSx_zstag(sig,i,j,k);}

    real dSy_dir_xstaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSy_xstag(sig,i,j,k);}
    real dSy_dir_ystaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSy_ystag(sig,i,j,k);}
    real dSy_dir_zstaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSy_zstag(sig,i,j,k);}

    real dSz_dir_xstaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSz_xstag(sig,i,j,k);}
    real dSz_dir_ystaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSz_ystag(sig,i,j,k);}
    real dSz_dir_zstaggered(const Sign sig, const int i, const int j, const int k) const
     {return dom->dSz_zstag(sig,i,j,k);}

    real dV_xstaggered(const int i, const int j, const int k) const {return dom->dV_xstag(i,j,k);}
    real dV_ystaggered(const int i, const int j, const int k) const {return dom->dV_ystag(i,j,k);}
    real dV_zstaggered(const int i, const int j, const int k) const {return dom->dV_zstag(i,j,k);}

    /* take care that all of the data bellow is passed in copy constructor */
    ScalarBoolP  vec;
    const    Domain * dom;

    const bool aliav;
};

#endif
