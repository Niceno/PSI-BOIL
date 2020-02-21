#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Field/ScalarInt/scalarint.h"
#include "../../Field/Vector/vector.h"
#include "../../Domain/domain.h"

////////////////
//            //
//  Topology  //
//            //
////////////////
/* A class for safe argument passing from interface tracking to other classes */
class Topology {
  public:
    Topology(Scalar * CLR, Scalar * NX, Scalar * NY, Scalar * NZ, 
             Scalar * ADENS, Vector * FS, ScalarInt * IFLAG) :
        clr(CLR),
        nx(NX),
        ny(NY),
        nz(NZ),
        fs(FS),
        adens(ADENS),
        iflag(IFLAG),
        stmp(*CLR->domain()),
        delta(*CLR->domain()),
        stmp2(*CLR->domain()) {
 
      stmp  = NX->shape();
      stmp2 = NX->shape();
      delta = NX->shape();

      mmax_ext = 100;
      tol_ext = 1e-7; 
    }
    
    ~Topology() {};

    void extrapolate(Scalar & sca, const Sign iext);
    void extrapolate(Scalar & sca, const Sign iext, const ScalarInt & eflag);
    
    inline int get_extrapolation_iters() const { return mmax_ext; }
    inline real get_extrapolation_tol() const  { return tol_ext; }
    inline void set_extrapolation_params(const int mnew, const real tolnew) {
      mmax_ext = mnew;
      tol_ext = tolnew;
      boil::oout<<"Topology::extrapolationparams: "<<mnew<<" "<<tolnew<<"\n";
    }

    Scalar * clr;
    Scalar * nx;
    Scalar * ny;
    Scalar * nz;
    Scalar * adens;
    Vector * fs;
    ScalarInt * iflag;

  private:
    int mmax_ext;
    real tol_ext;
    
    ScalarInt stmp;
    Scalar delta, stmp2;
};
#endif
