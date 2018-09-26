#ifndef PLOTGMV_H
#define PLOTGMV_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <climits>    /* INT_MAX, INT_MIN */

#include "../plot.h"
#include "../../Parallel/communicator.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Field/Vector/vector.h"

///////////////
//           //
//  PlotGMV  //
//           //
///////////////
class PlotGMV : public Plot {
  public:
    PlotGMV( const AsNodes asno=AsNodes::no(),
             const Buffers buff=Buffers::no() ) {
      sh=0;    if( buff == Buffers::yes() ) sh=1;
    }

    void plot(Domain &, const char *, const int);
    void plot(Body &, const char *, const int);
    void plot(const Scalar &, const char *, const int);
    void plot(const ScalarInt &, const char *, const int);
    void plot(const Vector &, const char *, const int);
    void plot(const Vector &, const Scalar &, const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &,
              const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const char *, const int);
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
    void plot_gmv_header(const char *, const int);
    void plot_gmv_domain(const Domain &);
    void plot_gmv_body  (const Body   &);
    void plot_gmv_start_scalars();
    void plot_gmv_end_scalars();
    void plot_gmv_scalar(const Domain &, const Scalar &, const char *);
    void plot_gmv_scalarint(const Domain &, const ScalarInt &, const char *);
    void plot_gmv_vector(const Domain &, const Vector &);
    void plot_gmv_footer();
    std::ofstream out;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: plot_gmv.h,v 1.20 2015/05/05 14:27:45 sato Exp $'/
+-----------------------------------------------------------------------------*/
