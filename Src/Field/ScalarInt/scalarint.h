#ifndef SCALARINT_H
#define SCALARINT_H

#include "../../Parallel/mpi_macros.h"
#include <sstream>
#include <fstream>
#include <cmath>
#include "../../Parallel/communicator.h"
#include "../../Domain/domain.h"
#include "../../Boundary/bndcnd.h"
#include "../Scalar/scalar_browsing.h"
#include "../../Global/global_malloc.h"
#include "../../Global/global_minmax.h"
#include "../../Global/global_name_file.h"
#include "../../LookUpTable/lookuptable.h"  
#include "../../Ravioli/sign.h"

#include "../Scalar/scalar_acc.h"
#include "../Scalar/scalar_aliases.h"

/////////////////
//             //
//  ScalarINT  //
//             //
/////////////////
class ScalarInt {
  public:
    /* global constructor */
    explicit ScalarInt(const Domain & d);
    explicit ScalarInt(const Domain & d, const char * n);
    explicit ScalarInt(const Domain & d, BndCnd & b);
    explicit ScalarInt(const ScalarInt & s); 
    explicit ScalarInt(const ScalarInt * s); // this creates an alias

    /* local constructors and functions */
    ScalarInt() : alias(true) {}; /* alias arbitrary here */

    void allocate(int ni, int nj, int nk); 

    /* destructor (which has been fixed) */
    ~ScalarInt();

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
	  
    int ** operator [] (const int i) const {return val[i];}

    int & operator () (int i, int j, int k) {
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
    int min() const;
    int max() const;
    int min_abs() const;
    int max_abs() const;
    int min_at() const;
    int max_at() const;

    void limit(const Range<real> & lim)
     {for_ijk(i,j,k) 
       {if(val[i][j][k] < lim.first()) val[i][j][k]=lim.first();
        if(val[i][j][k] > lim.last())  val[i][j][k]=lim.last();}}

    BndCnd & bc() const {return * bndcnd;}
    const Domain * domain() const {return dom;}
    
    void bnd_update();
    void bnd_update_nowall();
    void bnd_update_symmetry();

    void bnd_extract( const Dir d, int *** cp ) const;
    void bnd_insert ( const Dir d, int **  cp );

    /* the same as below but bndcnds are copied */
    ScalarInt & copy_shape(const Shape & a) {
      n_x=a.n_i;n_y=a.n_j;n_z=a.n_k;  o_x=a.o_i;o_y=a.o_j;o_z=a.o_k;
      s_x=a.s_i;s_y=a.s_j;s_z=a.s_k;  e_x=a.e_i;e_y=a.e_j;e_z=a.e_k;
      dom=a.dm;
      bc().replicate(*(a.bc));

      return *this;
    }

    /* operators */
    ScalarInt & operator = (const Shape & a)
     {n_x=a.n_i;n_y=a.n_j;n_z=a.n_k;  o_x=a.o_i;o_y=a.o_j;o_z=a.o_k;
      s_x=a.s_i;s_y=a.s_j;s_z=a.s_k;  e_x=a.e_i;e_y=a.e_j;e_z=a.e_k;
      bndcnd=a.bc; dom=a.dm;
      return *this;}

    /* mathematical operators */
    /* = */
    const ScalarInt & operator = (const ScalarInt & s) 
     {for_aijk(i,j,k) val[i][j][k]=s.val[i][j][k]; return *this;}
    const ScalarInt & operator = (const int & d) 
     {for_aijk(i,j,k) val[i][j][k]=d; return *this;}
    const Scalar & operator = (const char * c); 
    /* 
       += */
    const ScalarInt & operator += (const ScalarInt & s) 
     {for_aijk(i,j,k) val[i][j][k]+=s.val[i][j][k]; return *this;}
    const ScalarInt & operator += (const int & d)
     {for_aijk(i,j,k) val[i][j][k]+=d; return *this;}
    /* 
       -= */
    const ScalarInt & operator -= (const ScalarInt & s) 
     {for_aijk(i,j,k) val[i][j][k]-=s.val[i][j][k]; return *this;}
    const ScalarInt & operator -= (const int & d)
     {for_aijk(i,j,k) val[i][j][k]-=d; return *this;}
    /* 
       *= */
    const ScalarInt & operator *= (const ScalarInt & s) 
     {for_aijk(i,j,k) val[i][j][k]*=s.val[i][j][k]; return *this;}
    const ScalarInt & operator *= (const int & d) 
     {for_aijk(i,j,k) val[i][j][k]*=d; return *this;}
    /* 
       dot */
    real dot(const ScalarInt & s) 
     {real d=0.0; for_ijk(i,j,k) d+=val[i][j][k]*s.val[i][j][k]; 
      boil::cart.sum_real(&d); return d;}

    friend class VectorInt;

    const std::string & name() const {return nam;}

  protected:
    void deallocate();
 
    /* take care that all of the data bellow is passed in copy constructor */
    int ***     val;
    int         n_x, n_y, n_z;
    int         s_x, s_y, s_z;
    int         e_x, e_y, e_z;
    int         o_x, o_y, o_z;  // offset in x, y and z!
    const       Domain * dom;
    std::string nam;

    BndCnd * bndcnd;

    const bool alias;

  private:
    int miss; /* miss out of range variables in parallel version  */

};

#endif
