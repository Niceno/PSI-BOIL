#ifndef SCALAR_H
#define SCALAR_H

#include "../../Parallel/mpi_macros.h"
#include <sstream>
#include <fstream>
#include <cmath>
#include "../../Parallel/communicator.h"
#include "../../Domain/domain.h"
#include "../../Boundary/bndcnd.h"
#include "scalar_browsing.h"
#include "../../Global/global_malloc.h"
#include "../../Global/global_minmax.h"
#include "../../Global/global_name_file.h"
#include "../../LookUpTable/lookuptable.h"  
#include "../../Ravioli/sign.h"

#include "scalar_acc.h"
#include "scalar_aliases.h"

//////////////
//          //
//  Scalar  //
//          //
//////////////
class Scalar {
  public:
    /* global constructor */
    explicit Scalar(const Domain & d);
    explicit Scalar(const Domain & d, const char * n);
    explicit Scalar(const Domain & d, BndCnd & b);
    explicit Scalar(const Scalar & s); 
    explicit Scalar(const Scalar * s); // this creates an alias

    /* local constructors and functions */
    Scalar() : alias(true) {}; /* alias arbitrary here */

    void allocate(int ni, int nj, int nk); 

    /* destructor (which has been fixed) */
    ~Scalar();

    /* used for (re)dimensioning arrays */
    /* (for helping arrays in linear solvers) */
    const Shape shape() const
     {return Shape(n_x,n_y,n_z,  o_x,o_y,o_z,  s_x,s_y,s_z,  e_x,e_y,e_z,
                   bndcnd, dom);}

    int ni() const {return n_x;}      
    int nj() const {return n_y;}      
    int nk() const {return n_z;}      
    int ox() const {return o_x;}      
    int oy() const {return o_y;}      
    int oz() const {return o_z;}      
    int si() const {return s_x;}      
    int sj() const {return s_y;}      
    int sk() const {return s_z;}      
    int ei() const {return e_x;}      
    int ej() const {return e_y;}      
    int ek() const {return e_z;}      

    /* change (set) offsets */
    void ox(const int o) {o_x=o;}      
    void oy(const int o) {o_y=o;}      
    void oz(const int o) {o_z=o;}
      
    /* cell coordinates (could these be pointers?) */
    real dxc(const int i) const {return dom->dxc(i);}
    real dyc(const int j) const {return dom->dyc(j);}
    real dzc(const int k) const {return dom->dzc(k);}

    real dxw(const int i) const {return dom->dxw(i);}
    real dxe(const int i) const {return dom->dxe(i);}
    real dys(const int j) const {return dom->dys(j);}
    real dyn(const int j) const {return dom->dyn(j);}
    real dzb(const int k) const {return dom->dzb(k);}
    real dzt(const int k) const {return dom->dzt(k);}

    /* cell centre coordinates (needed for Formula evaluation */
    real xc(const int i) const {return dom->xc(i);}
    real yc(const int j) const {return dom->yc(j);}
    real zc(const int k) const {return dom->zc(k);}

    /* node coordinates */
    real xn(const int i) const {return dom->xn(i);}
    real yn(const int j) const {return dom->yn(j);}
    real zn(const int k) const {return dom->zn(k);}

    /* carefull: these return local logical coordinates (for Body) */
    int i(const real x) const;
    int j(const real y) const;
    int k(const real z) const;
    int i(const real x, const real t) const;
    int j(const real y, const real t) const;
    int k(const real z, const real t) const;
    int im(const real x, const real t) const;
    int jm(const real y, const real t) const;
    int km(const real z, const real t) const;
    int ip(const real x, const real t) const;
    int jp(const real y, const real t) const;
    int kp(const real z, const real t) const;
    int aim(const real x, const real t) const;
    int ajm(const real y, const real t) const;
    int akm(const real z, const real t) const;
    int aip(const real x, const real t) const;
    int ajp(const real y, const real t) const;
    int akp(const real z, const real t) const;

    real dSx(const int i, const int j, const int k) const 
     {return dom->dSx(i,j,k);}
    real dSy(const int i, const int j, const int k) const 
     {return dom->dSy(i,j,k);}
    real dSz(const int i, const int j, const int k) const 
     {return dom->dSz(i,j,k);}

    real dSx(const Sign sig, const int i, const int j, const int k) const 
     {return dom->dSx(sig,i,j,k);}
    real dSy(const Sign sig, const int i, const int j, const int k) const 
     {return dom->dSy(sig,i,j,k);}
    real dSz(const Sign sig, const int i, const int j, const int k) const 
     {return dom->dSz(sig,i,j,k);}

    /* cell volume */
    real dV(const int i, const int j, const int k) const 
     {return dom->dV(i,j,k);}

    void exchange        (const int dir = -1) const; // should it be const???
    void exchange_avg    (const int dir)      const; // should it be const???
    void exchange_all    (const int dir = -1) const; // should it be const???
    void exchange_all_avg(const int dir = -1) const; // should it be const???
    void exchange        (const int * i, const int dir = -1) const;

    void one2all     (const int one);
	  
    real ** operator [] (const int i) const {return val[i];}

    real & operator () (int i, int j, int k) {
      if( dom->contains_IJK(i,j,k) ) {
        dom->locals(&i,&j,&k); // change i,j,k
        return val[i][j][k];   // with changed i,j,k
      }
      else
        return miss;
    }
    
    void save(const char *, const int);
    void load(const char *, const int);
    void rm  (const char *, const int);
    void save(std::ofstream &);
    void load(std::ifstream &);

    /* min, max, ... */
    real min() const;
    real max() const;
    real min_abs() const;
    real max_abs() const;
    real min_at() const;
    real max_at() const;

    /* capital letters mean global logical coordinates */
    real average_I(int I) const; 
    real average_J(int J) const; 
    real average_K(int K) const; 
    real average_IJ(int I, int J) const; 
    real average_IK(int I, int K) const; 
    real average_JK(int J, int K) const; 

    /* lower-case letters mean local logical coordinates */
    real average_i(int i) const; 
    real average_j(int j) const; 
    real average_k(int k) const; 
    real average_ij(int i, int j) const; 
    real average_ik(int i, int k) const; 
    real average_jk(int j, int k) const; 

    void limit(const Range<real> & lim)
     {for_ijk(i,j,k) 
       {if(val[i][j][k] < lim.first()) val[i][j][k]=lim.first();
        if(val[i][j][k] > lim.last())  val[i][j][k]=lim.last();}}

    real integral() const
     {real tot=0.0;
      for_ijk(i,j,k) tot += val[i][j][k] * dV(i,j,k);
      boil::cart.sum_real(&tot); 
      return tot;}

    BndCnd & bc() const {return * bndcnd;}
    const Domain * domain() const {return dom;}
    
    void grad( const int i, const int j, const int k,
               real * nx, real * ny, real * nz) const;
    void grad_abs( const int i, const int j, const int k,
                   real * nx, real * ny, real * nz) const;

    void bnd_update();
    void bnd_update_nowall();
    void bnd_update_symmetry();
    void bnd_grad_update(const Comp &);

    void bnd_extract( const Dir d, real *** cp ) const;
    void bnd_insert ( const Dir d, real **  cp );

    /* look up for tabulated values */
    void look_up(const Scalar      & sca,
                 const LookUpTable & tab,
                 const Column      & col0, 
                 const Column      & col1) {
      for_avijk(sca, i, j, k) 
        val[i][j][k] = tab.look_up(sca[i][j][k], col0, col1);
     }

    /* the same as below but bndcnds are copied */
    Scalar & copy_shape(const Shape & a) {
      n_x=a.n_i;n_y=a.n_j;n_z=a.n_k;  o_x=a.o_i;o_y=a.o_j;o_z=a.o_k;
      s_x=a.s_i;s_y=a.s_j;s_z=a.s_k;  e_x=a.e_i;e_y=a.e_j;e_z=a.e_k;
      dom=a.dm;
      bc().replicate(*(a.bc));

      return *this;
    }

    /* operators */
    Scalar & operator = (const Shape & a)
     {n_x=a.n_i;n_y=a.n_j;n_z=a.n_k;  o_x=a.o_i;o_y=a.o_j;o_z=a.o_k;
      s_x=a.s_i;s_y=a.s_j;s_z=a.s_k;  e_x=a.e_i;e_y=a.e_j;e_z=a.e_k;
      bndcnd=a.bc; dom=a.dm;
      return *this;}

    /* mathematical operators */
    /* = */
    const Scalar & operator = (const Scalar & s) 
     {for_aijk(i,j,k) val[i][j][k]=s.val[i][j][k]; return *this;}

    const Scalar & operator = (const real & d) 
     {for_aijk(i,j,k) val[i][j][k]=d; return *this;}
    const Scalar & operator = (const char * c); 
    /* 
       += */
    const Scalar & operator += (const Scalar & s) 
     {for_aijk(i,j,k) val[i][j][k]+=s.val[i][j][k]; return *this;}
    const Scalar & operator += (const real & d)
     {for_aijk(i,j,k) val[i][j][k]+=d; return *this;}
    /* 
       -= */
    const Scalar & operator -= (const Scalar & s) 
     {for_aijk(i,j,k) val[i][j][k]-=s.val[i][j][k]; return *this;}
    const Scalar & operator -= (const real & d)
     {for_aijk(i,j,k) val[i][j][k]-=d; return *this;}
    /* 
       *= */
    const Scalar & operator *= (const Scalar & s) 
     {for_aijk(i,j,k) val[i][j][k]*=s.val[i][j][k]; return *this;}
    const Scalar & operator *= (const real & d) 
     {for_aijk(i,j,k) val[i][j][k]*=d; return *this;}
    /* 
       /= */
    const Scalar & operator /= (const Scalar & s) 
     {for_aijk(i,j,k) val[i][j][k]/=s.val[i][j][k]; return *this;}
    const Scalar & operator /= (const real & d) 
     {for_aijk(i,j,k) val[i][j][k]/=d; return *this;}
    /* 
       dot */
    real dot(const Scalar & s) 
     {real d=0.0; for_ijk(i,j,k) d+=val[i][j][k]*s.val[i][j][k]; 
      boil::cart.sum_real(&d); return d;}
    /*
       sum */
    real sum() {
      real r(0.);
      for_ijk(i,j,k)
        r += val[i][j][k];
      boil::cart.sum_real(&r);
      return r;
    }
    /*
       sum */
    real sum_abs() {
      real r(0.);
      for_ijk(i,j,k)
        r += fabs(val[i][j][k]);
      boil::cart.sum_real(&r);
      return r;
    }

    /*-------------------------------------+
    |  accelerated mathematical operators  |
    +-------------------------------------*/
    /* 
       = alfa x */
    const Scalar & operator = (const alfa_x & c)
     {for_aijk(i,j,k) val[i][j][k] = c.alfa * c.x[i][j][k]; return *this;}
    /* 
       = x alfa */
    const Scalar & operator = (const x_alfa & c)
     {for_aijk(i,j,k) val[i][j][k] = c.x[i][j][k] * c.alfa; return *this;}
    /* 
       = x / alfa */
    const Scalar & operator = (const x_d_alfa & c)
     {for_aijk(i,j,k) val[i][j][k] = c.x[i][j][k] / c.alfa; return *this;}
    /* 
       += alfa x */
    const Scalar & operator += (const alfa_x & c)
     {for_aijk(i,j,k) val[i][j][k] += c.alfa * c.x[i][j][k]; return *this;}
    /* 
       -= alfa x */
    const Scalar & operator -= (const alfa_x & c)
     {for_aijk(i,j,k) val[i][j][k] -= c.alfa * c.x[i][j][k]; return *this;}
    /* 
       = y + alfa x */
    const Scalar & operator = (const y_p_alfa_x & c)
     {for_aijk(i,j,k) val[i][j][k] = c.y[i][j][k] + c.c.alfa * c.c.x[i][j][k]; 
      return *this;}
    /* 
       = y - alfa x */
    const Scalar & operator = (const y_m_alfa_x & c)
     {for_aijk(i,j,k) val[i][j][k] = c.y[i][j][k] - c.c.alfa * c.c.x[i][j][k]; 
      return *this;}
    /* 
       = y - A x */
    const Scalar & operator = (const y_m_A_x & c);
    /* 
       = A x */
    const Scalar & operator = (const A_x & c);
    /* 
       = x + y */
    const Scalar & operator = (const x_p_y & c)
     {for_aijk(i,j,k) val[i][j][k] = c.x[i][j][k] + c.y[i][j][k]; return *this;}
    /* 
       = x - y */
    const Scalar & operator = (const x_m_y & c)
     {for_aijk(i,j,k) val[i][j][k] = c.x[i][j][k] - c.y[i][j][k]; return *this;}
    /* 
       = x * y (NOT a dot!) */
    const Scalar & operator = (const x_y & c)
     {for_aijk(i,j,k) val[i][j][k] = c.x[i][j][k] * c.y[i][j][k]; return *this;}
    /* 
       = x / y */
    const Scalar & operator = (const x_d_y & c)
     {for_aijk(i,j,k) val[i][j][k] = c.x[i][j][k] / c.y[i][j][k]; return *this;}
    /* 
       += x * y (NOT a dot!) */
    const Scalar & operator += (const x_y & c)
     {for_aijk(i,j,k) val[i][j][k] += c.x[i][j][k] * c.y[i][j][k]; return *this;}
    /* 
       += x / y */
    const Scalar & operator += (const x_d_y & c)
     {for_aijk(i,j,k) val[i][j][k] += c.x[i][j][k] / c.y[i][j][k]; return *this;}
    /* 
       -= x * y (NOT a dot!) */
    const Scalar & operator -= (const x_y & c)
     {for_aijk(i,j,k) val[i][j][k] -= c.x[i][j][k] * c.y[i][j][k]; return *this;}
    /* 
       -= x / y */
    const Scalar & operator -= (const x_d_y & c)
     {for_aijk(i,j,k) val[i][j][k] -= c.x[i][j][k] / c.y[i][j][k]; return *this;}

    friend class Vector;

    const std::string & name() const {return nam;}

  protected:
    void deallocate();
 
    /* take care that all of the data bellow is passed in copy constructor */
    real ***    val;
    int         n_x, n_y, n_z;
    int         s_x, s_y, s_z;
    int         e_x, e_y, e_z;
    int         o_x, o_y, o_z;  // offset in x, y and z!
    const       Domain * dom;
    std::string nam;

    BndCnd * bndcnd;

    const bool alias;

  private:
    real miss; /* miss out of range variables in parallel version  */

};

#endif
