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

    // unused functions (still necessary because of pure virtual function)
    void plot(const char *, const int, Times * t,
                      const Scalar * s1,
                      const Scalar * s2 = NULL,
                      const Scalar * s3 = NULL,
                      const Scalar * s4 = NULL,
                      const Scalar * s5 = NULL,
                      const Scalar * s6 = NULL,
                      const Scalar * s7 = NULL,
                      const Scalar * s8 = NULL,
                      const Scalar * s9 = NULL) {
       boil::oout<<"Error! These arguments are for TECMPI\n"; exit(0);};
     void plot(const char *, const int, Times * t, const Vector *,
                      const Scalar * s1 = NULL,
                      const Scalar * s2 = NULL,
                      const Scalar * s3 = NULL,
                      const Scalar * s4 = NULL,
                      const Scalar * s5 = NULL,
                      const Scalar * s6 = NULL,
                      const Scalar * s7 = NULL,
                      const Scalar * s8 = NULL,
                      const Scalar * s9 = NULL) {
       boil::oout<<"Error! These arguments are for TECMPI\n"; exit(0);};
     void read(const char *, const int, Times * t, Vector *,
                      Scalar * s1 = NULL,
                      Scalar * s2 = NULL,
                      Scalar * s3 = NULL,
                      Scalar * s4 = NULL,
                      Scalar * s5 = NULL,
                      Scalar * s6 = NULL,
                      Scalar * s7 = NULL,
                      Scalar * s8 = NULL,
                      Scalar * s9 = NULL) {
       boil::oout<<"Error! This function is valid only  for TECMPI\n"; exit(0);};
     void read(const char *, const int, Times * t,
                      Scalar * s1 = NULL,
                      Scalar * s2 = NULL,
                      Scalar * s3 = NULL,
                      Scalar * s4 = NULL,
                      Scalar * s5 = NULL,
                      Scalar * s6 = NULL,
                      Scalar * s7 = NULL,
                      Scalar * s8 = NULL,
                      Scalar * s9 = NULL) {
       boil::oout<<"Error! This function is valid only  for TECMPI\n"; exit(0);};

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
