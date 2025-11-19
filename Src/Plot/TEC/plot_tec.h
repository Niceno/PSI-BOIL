#ifndef PLOTTEC_H
#define PLOTTEC_H

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
//  PlotTEC  //
//           //
///////////////
class PlotTEC : public Plot {
  public:
    PlotTEC( const AsNodes asno=AsNodes::no(),
             const Buffers buff=Buffers::no() ) {
      nodal=0; if( asno == AsNodes::yes() ) nodal=1;
      sh=0;    if( buff == Buffers::yes() ) sh=1;
      b_plot_body=true;
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
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &,
              const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int);
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const char *, const int);
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
    void set_plot_body(bool b){
      b_plot_body = b;
      boil::oout<<"PlotTEC:set_plot_body= "<<b<<"\n";
    }

  private:
    void plot_tec_header   (const Domain &, const char *, const int);
    void plot_tec_header   (const char   *, const int);
    void plot_tec_prologue (const std::vector<std::string> & vnames);
    void plot_tec_domain   (const Domain &);
    void plot_tec_body     (const Body   &, std::vector<int> &);
    void plot_tec_scalar   (const Domain &, const Scalar &);
    void plot_tec_scalarint(const Domain &, const ScalarInt &);
    void plot_tec_vector   (const Domain &, const Vector &);
    void plot_tec_footer   ();
    std::ofstream out;
    bool b_plot_body;
};

#endif
