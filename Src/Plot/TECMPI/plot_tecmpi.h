#ifndef PLOTTECMPI_H
#define PLOTTECMPI_H
#define TECIOMPI

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <climits>    /* INT_MAX, INT_MIN */
#include <vector>
#include <memory>

#include "../plot.h"
#include "../../Parallel/communicator.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Field/ScalarInt/scalarint.h"
#include "../../Field/Vector/vector.h"
#include "../../Global/global_buffer_width.h"
#include "teciompisrc/TECIO.h"

#ifndef NULL
#define NULL 0
#endif

#define VAR_TYPE real
#define TECIOMPI

using namespace std;

///////////////
//           //
//  PlotTEC  //
//           //
///////////////
class PlotTECMPI : public Plot {
  public:
    PlotTECMPI(const AsNodes asno=AsNodes::no(),
          const Buffers buff=Buffers::no()) {

      nodal=0; if( asno == AsNodes::yes() ) nodal=1;
      sh=0;    if( buff == Buffers::yes() ) sh=1;
      b_plot_body=false;
      Debug      = 1;
      isDouble  = 0;
      FileType   = 0;
      fileFormat = 1; // SZPLT; .PLT not supported for partitioned zones
      I          = 0; // Used to track return codes
      mainRank   = 0;
      commSize   = boil::cart.nproc();
      commRank   = boil::cart.iam();
      mpiComm    = boil::cart.world();
      numZones   = commSize;
      BW         = boil::BW;
    }
    void plot(const char *, const int, Times * t,
                      const Scalar * s1,
                      const Scalar * s2 = NULL,
                      const Scalar * s3 = NULL,
                      const Scalar * s4 = NULL,
                      const Scalar * s5 = NULL,
                      const Scalar * s6 = NULL,
                      const Scalar * s7 = NULL,
                      const Scalar * s8 = NULL,
                      const Scalar * s9 = NULL);
    void plot(const char *, const int, Times * t, const Vector *,
                      const Scalar * s1 = NULL,
                      const Scalar * s2 = NULL,
                      const Scalar * s3 = NULL,
                      const Scalar * s4 = NULL,
                      const Scalar * s5 = NULL,
                      const Scalar * s6 = NULL,
                      const Scalar * s7 = NULL,
                      const Scalar * s8 = NULL,
                      const Scalar * s9 = NULL);
    void plot(Domain &, const char *, const int, Times * t = NULL);
    void plot(const Pathline &, const char *, const int, Times * t = NULL);
    void read(const char *, const int, Times * t, Vector *,
                      Scalar * s1 = NULL,
                      Scalar * s2 = NULL,
                      Scalar * s3 = NULL,
                      Scalar * s4 = NULL,
                      Scalar * s5 = NULL,
                      Scalar * s6 = NULL,
                      Scalar * s7 = NULL,
                      Scalar * s8 = NULL,
                      Scalar * s9 = NULL);

    // unused functions (still necessary because of pure virtual function)
    void plot(Body &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const ScalarInt &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const ScalarInt &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void plot(const Vector &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const Scalar &, const Scalar &,
              const Scalar &, const Scalar &, const char *, const int, Times * t = NULL) {
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};
    void set_plot_body(bool b){
         boil::oout<<"Error! Check arguments for TECMPI\n"; exit(0);};

  private:
    void plot_tecmpi_set_domain(const Domain * d);
    string plot_tecmpi_fname(const char * nam, const int i);
    string plot_tecmpi_variables(const vector<string>& v);
    void plot_tecmpi_write (const Domain * d, int32_t & zone, int & i);
    void plot_tecmpi_write (const Vector * v, int32_t & zone, int & i);
    void plot_tecmpi_write (const Scalar * s, int32_t & zone, int & i);
    void plot_tecmpi_tecOpenInit(const string & s1, const string & s2);
    void plot_tecmpi_tecCreateMap(vector<int32_t> & v1, vector<int32_t> & v2, 
           vector<int32_t> & v3, int32_t & i1, int32_t & i2, const int &i,
           const Times * t = NULL);
    void copy_valCell(const std::unique_ptr<float[]>& values, Scalar *s);
    void copy_valCell(const std::unique_ptr<float[]>& values,
                              Vector *v, const Comp &m);

    ofstream out;
    bool b_plot_body;

    // for teciompi
    INTEGER4 Debug,isDouble,FileType,fileFormat,I;
    INTEGER4 mainRank;
    par_comm mpiComm;
    int      commSize, commRank;
    int      numZones, numVars;
    int      XDIM, YDIM, ZDIM, XDIM_C, YDIM_C, ZDIM_C;
    int32_t res;
    void*   fileHandle;
    int BW;

};

#endif
