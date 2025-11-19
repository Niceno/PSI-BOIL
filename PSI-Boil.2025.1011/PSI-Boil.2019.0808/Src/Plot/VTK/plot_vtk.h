#ifndef PLOTVTK_H
#define PLOTVTK_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <climits>    /* INT_MAX, INT_MIN */

#include "../plot.h"
#include "../../Parallel/communicator.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Field/ScalarInt/scalarint.h"
#include "../../Field/Vector/vector.h"

///////////////
//           //
//  PlotVTK  //
//           //
///////////////
class PlotVTK : public Plot {
  public:
    PlotVTK( const AsNodes asno=AsNodes::no(),
             const Buffers buff=Buffers::no() ) {
      sh=0;    if( buff == Buffers::yes() ) sh=1;
    }

    void plot(Domain &, const char *, const int);
    void plot(Body &, const char *, const int = -1);
    void plot(const Scalar &, const char *, const int);
    void plot(const ScalarInt &, const char *, const int);
    void plot(const Vector &, const char *, const int);
    void plot(const Vector &, const Scalar &, const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &, 
              const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int); 
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int); 
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const char *, const int); 
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &,
              const char *, const int); 
    void plot(const Scalar &, const Scalar &, 
              const char *, const int);
    void plot(const Scalar &, const ScalarInt &, 
              const char *, const int);
    void plot(const Scalar &, const Scalar &, const Scalar &,
              const char *, const int);
    void plot(const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const char *, const int);
    void plot(const Scalar &, const Scalar &, const Scalar &, 
              const Scalar &, const Scalar &, 
              const char *, const int);
    void plot(const Scalar &, const Scalar &, const Scalar &, 
              const Scalar &, const Scalar &, const Scalar &,
              const char *, const int);

  private:
    void plot_vtk_header   (const Domain &, const char *, const int);
    void plot_vtk_domain   (const Domain &);
    void plot_vtk_scalar   (const Domain &, const Scalar &, const char *);
    void plot_vtk_scalarint(const Domain &, const ScalarInt &, const char *);
    void plot_vtk_vector   (const Domain &, const Vector &);
    void plot_vtk_footer   ();
    std::ofstream out;
};

#endif
