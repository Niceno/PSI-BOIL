#include "marching_squares.h"

MarchingSquares::MarchingSquares(const Comp & M, const Scalar * CLR,
                                 Scalar * PHI, Scalar * ADENS,
                                 const real CLRSURF) :
  Heaviside(CLR,PHI,ADENS,CLRSURF),
  stmp(*(CLR->domain())),
  perpendicular(M)

 {

  stmp = CLR->shape();

  /* positional offset */
  ofx = -1;
  ofy = -1;
  ofz = -1;

  /* nodal values */
  nodalvals = CLR->shape();
  int ofx_alloc(1), ofy_alloc(1), ofz_alloc(1);
  if       (perpendicular == Comp::i()) {
    ofx_alloc--;
    ofx++;
  } else if(perpendicular == Comp::j()) {
    ofy_alloc--;
    ofy++;
  } else if(perpendicular == Comp::k()) {
    ofz_alloc--;
    ofz++;
  } else {
    boil::oout<<"Marching Squares direction not properly set! Exiting."
              <<boil::endl;
    exit(0);
  }
  nodalvals.allocate(CLR->ni()+ofx_alloc, 
                     CLR->nj()+ofy_alloc,
                     CLR->nk()+ofz_alloc);
  nodalvals.ox(ofx_alloc);
  nodalvals.oy(ofy_alloc);
  nodalvals.oz(ofz_alloc);

  return;
}  
