#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
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
    Topology(Scalar * NX, Scalar * NY, Scalar * NZ, 
             Scalar * ADENS, Vector * FS, ScalarInt * IFLAG) :
        nx(NX),
        ny(NY),
        nz(NZ),
        fs(FS),
        adens(ADENS),
        iflag(IFLAG) { }
    ~Topology() {};

    Scalar * nx;
    Scalar * ny;
    Scalar * nz;
    Scalar * adens;
    Vector * fs;
    ScalarInt * iflag;
};
#endif
