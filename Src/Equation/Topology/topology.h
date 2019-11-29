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
        iflag(IFLAG) { }
    ~Topology() {};

    Scalar * clr;
    Scalar * nx;
    Scalar * ny;
    Scalar * nz;
    Scalar * adens;
    Vector * fs;
    ScalarInt * iflag;
};
#endif
