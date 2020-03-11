#include "linear.h"

/***************************************************************************//**
*  \brief Parent class class constructor for all linear algebraic solvers.
*
*  Stores the domain on which the solver is created and sets the 
*  Preconditioner. The defualt Preconditioner is "di", diagonal.
*
*  \note
*    It does not allocate the memory for children solvers, since they are
*    different for each child. 
*******************************************************************************/
Linear::Linear(const Domain & d, const Prec & pc) 
 : /* p (), 
   q (), 
   r (), 
   s (), 
   u (), 
   z (), 
   p_(), 
   q_(), 
   r_(), 
   u_(), 
   v_(), */ 
   dom(& d) {

  assert(d.ni() > 0);
  assert(d.nj() > 0);
  assert(d.nk() > 0);

  if(pc == Prec::di ()) prec = new Diagonal(d);
  if(pc == Prec::ic0()) prec = new IncompleteCholesky0(d);
  if(pc == Prec::ic2()) prec = new IncompleteCholesky2(d);
  if(pc == Prec::ic3()) prec = new IncompleteCholesky3(d);
}
