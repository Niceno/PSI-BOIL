#ifndef PLOT_H
#define PLOT_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include "../Parallel/Out/out.h"
#include "../Timer/timer.h"
#include "../Global/global_name_file.h"
#include "../Body/body.h"
#include "../Ravioli/buffers.h"
#include "../Ravioli/asnodes.h"
#include "../Ravioli/box.h"
#include "../Domain/domain.h"
#include "../Field/ScalarInt/scalarint.h"

class Scalar;
class Vector;
class ScalrInt;

////////////
//        //
//  Plot  //
//        //
////////////
class Plot {
  public:
    Plot( const AsNodes asno=AsNodes::no(),
          const Buffers buff=Buffers::no() ) {} 

    virtual void plot(Domain &, const char *, const int = -1) = 0;
    virtual void plot(Body &, const char *, const int = -1) = 0;
    virtual void plot(const Scalar &, const char *, const int = -1) = 0;
    virtual void plot(const ScalarInt &, const char *, const int = -1) = 0;
    virtual void plot(const Vector &, const char *, const int = -1) = 0;
    virtual void plot(const Vector &, const Scalar &, const char *, const int) 
                 = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &,
                      const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const Scalar &, const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, 
                      const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const char *, const int) = 0;
    virtual void plot(const Vector &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const Scalar &, const Scalar &, const Scalar &,
                      const Scalar &,
                      const char *, const int) = 0;
    virtual void plot(const Scalar &, const Scalar &, 
                      const char *, const int = -1) = 0;
    virtual void plot(const Scalar &, const ScalarInt &, 
                      const char *, const int = -1) = 0;
    virtual void plot(const Scalar &, const Scalar &, const Scalar &,
                      const char *, const int = -1) = 0;
    virtual void plot(const Scalar &, const Scalar &, 
                      const Scalar &, const Scalar &,
                      const char *, const int = -1) = 0;
    virtual void plot(const Scalar &, const Scalar &, const Scalar &, 
                      const Scalar &, const Scalar &, 
                      const char *, const int) = 0;
    virtual void plot(const Scalar &, const Scalar &, const Scalar &, 
                      const Scalar &, const Scalar &, const Scalar &,
                      const char *, const int) = 0;
    virtual void set_plot_body(bool b) = 0;
  protected:
    const Domain * dom;
    int   sh;           /* shift */
    int   nodal;        /* if 1 plot as nodes */
};	

/*-----------------+
|  global plotter  |
+-----------------*/
namespace boil {
  extern Plot * plot;
}

#endif
