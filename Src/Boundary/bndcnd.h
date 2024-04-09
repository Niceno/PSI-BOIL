#ifndef BNDCND_H
#define BNDCND_H

#include "../Parallel/mpi_macros.h"
#include <vector>
#include "../Ravioli/dir.h"
#include "../Ravioli/comp.h"
#include "../Ravioli/range.h"
#include "../Global/global_precision.h"
#include "../Domain/domain.h"
#include "../Formula/formula.h"
#include "../Body/body.h"
#include "../Global/global_endl.h"

/***************************************************************************//**
*  \brief Ravioli class for safer argument passing in BndCnd.
*
*  Used for setting the right boundary condition type in BndCnd and Variable. 
*  The names of the member functions are self-explanatory - they closely follow 
*  the boundary condition types one can use in a simulation.
*
*  For example, if we wanted to add Neumann boundary condition in "kmin"
*  direction of variable "t", we would us the code:
*  \code
*    Scalar t(d);
*    ...
*    t.bc()->add( BndCnd( Dir::kmin(), BndType::neumann() ) );
*  \endcode
*******************************************************************************/

///////////////
//           //
//  BndType  //
//           //
///////////////
class BndType {
  public:
    static const BndType undefined() {return BndType(0 );} 
    static const BndType dirichlet() {return BndType(1 );} // p t
    static const BndType neumann()   {return BndType(2 );} // p t
    static const BndType periodic()  {return BndType(3 );} // m p t
    static const BndType inlet()     {return BndType(4 );} // m t
    static const BndType outlet()    {return BndType(5 );} // m t
    static const BndType wall()      {return BndType(6 );} // m
    static const BndType symmetry()  {return BndType(7 );} // m p t
    static const BndType insert()    {return BndType(8 );} // m p t
    static const BndType convective(){return BndType(9 );} // m p t
    static const BndType pseudo()    {return BndType(10);} // m p t
    /* check << operator if above are changing */
    
    BndType(const BndType & bc) {val = bc.val;}

    //! Comapers if boundary types are equal.
    bool operator == (const BndType & other) const 
      {return val == other.val;}

    //! Comapers if boundary types are different.
    bool operator != (const BndType & other) const 
      {return val != other.val;}

    //! Prints the boundary value name.
    friend std::ostream & 
      operator << (std::ostream & os, const BndType & bt);

  private:
    /* prevent creation of new b.c. */
    explicit BndType(const int i) {val = i;}
    
    int val;
};	

/***************************************************************************//**
*  \brief Ravioli class for safer specification of arguments representing
*         boundary condition values.
*
*  It may hold either a real value if boundary condition is specified with a 
*  single constant, or with a string (character array), if boundary condition
*  is specified with an analytical expression.
*******************************************************************************/

//////////////
//          //
//  BndVal  //
//          //
//////////////
struct BndVal {
  BndVal()          : r(0.0), c(NULL) {}
  BndVal(real   r_) : r(r_),  c(NULL) {}
  BndVal(char * c_) : r(0.0), c(c_) {}
  real   r;
  char * c;   
};

//////////////
//  BndCnd  //
//////////////
class BndCnd {
  public:
    BndCnd( const Domain & d );  
    BndCnd( const Dir & dr, const BndType & bt, const BndVal & u = 0.0, 
                                                const BndVal & v = 0.0, 
                                                const BndVal & w = 0.0); 
    BndCnd( const Dir        & dr,  
            const Range<int> & j_r,
            const Range<int> & k_r, const BndType & bt, const BndVal & u = 0.0,
                                                        const BndVal & v = 0.0, 
                                                        const BndVal & w = 0.0 ); 
    BndCnd( const Range<int> & i_r, 
            const Dir        & dr, 
            const Range<int> & k_r, const BndType & bt, const BndVal & u = 0.0,
                                                        const BndVal & v = 0.0, 
                                                        const BndVal & w = 0.0 ); 
    BndCnd( const Range<int> & i_r, 
            const Range<int> & j_r,
            const Dir        & dr,  const BndType & bt, const BndVal & u = 0.0,
                                                        const BndVal & v = 0.0, 
                                                        const BndVal & w = 0.0 ); 
    BndCnd( const Range<int> & i_r, 
            const Range<int> & j_r,
            const Range<int> & k_r,
            const Dir        & dr,  const BndType & bt, const BndVal & u = 0.0,
                                                        const BndVal & v = 0.0, 
                                                        const BndVal & w = 0.0 ); 

    int si() const {return ir.first();}      
    int sj() const {return jr.first();}      
    int sk() const {return kr.first();}      
    int ei() const {return ir.last();}      
    int ej() const {return jr.last();}      
    int ek() const {return kr.last();}      

    void si(const int i) {ir.first(i);}      
    void sj(const int j) {jr.first(j);}      
    void sk(const int k) {kr.first(k);}      
    void ei(const int i) {ir.last(i);}      
    void ej(const int j) {jr.last(j);}      
    void ek(const int k) {kr.last(k);}      

    bool contains_ijk(const int i, const int j, const int k) {
      return ir.contains(i) && jr.contains(j) && kr.contains(k);
    }

    /* explicit operator was necessary for Intel C++ compiler */
    const BndCnd & operator = (const BndCnd & o)
     {typ    = o.typ;
      val[0] = o.val[0]; val[1] = o.val[1]; val[2] = o.val[2];
      return * this;}

    const Domain * domain() const {return dom;}

    /* manipulate with boundary conditions the same as in Vector */
          BndCnd &   at     (const int b) {return section[b];}
    void             add    (BndCnd bc);
    void             modify (BndCnd bc);
          BndType &  type   (const int b); 
    char *           formula(const int b, const Comp m=Comp::u()) const;
    real             value  (const int b, const Comp m=Comp::u()) const;
    real             value  (const int b, 
                             const real x, const real y, const real z,
                             const Comp m=Comp::u()) const;
    int              count() const;
    const Dir &      direction  (const int b) const;
    bool             type       (const Dir & dir, const BndType & bt) const;
    bool             type_decomp(const int b) const;
    bool             type_decomp(const Dir dir) const;
    bool             type_here  (const Dir & dir, const BndType & bt) const;

    bool exists(const int b) const;
    bool exists(const BndType & bt) const;  
    int  index (const Dir & dir, const BndType & bt) const;

    friend std::ostream & 
      operator << (std::ostream & os, const BndCnd & bc);

  private:
    BndType        typ;
    Dir            dir;
    Range<int>     ir, jr, kr; /* i,j,k range */
    BndVal         val[3];
    const Domain * dom; /* to check boundary condition extents */

    std::vector<BndCnd> section;
};	

#endif
