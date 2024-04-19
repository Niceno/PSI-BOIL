#ifndef EQUATION_H
#define EQUATION_H

#include "../Parallel/mpi_macros.h"
#include "../Field/Vector/vector.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/ScalarInt/scalarint.h"
#include "../Global/global_limiter.h"
#include "../Ravioli/timescheme.h"
#include "../SimulationTime/simulation_time.h"
#include "../Domain/domain.h"
#include "../Solver/Linear/Krylov/krylov.h"
#include "../Matter/matter.h"

////////////////
//            //
//  Equation  //
//            //
////////////////
class Equation {
  public:
    Equation(const Domain * d, 
             const Times  * t, 
             Matter * f, 
             Matter * s, 
             Linear * sm) : 
      dom(d), time(t), flu(f), sol(s), solver(sm),
      conv_ts(TimeScheme::adams_bashforth()),
      diff_ts(TimeScheme::crank_nicolson()),
      lim(ConvScheme::superbee()) {}  

    virtual void discretize(const Scalar * diff_eddy = NULL) = 0;

    virtual void diffusion_set (const TimeScheme & ts) {diff_ts = ts; discretize();}
    void convection_set(const TimeScheme & ts) {conv_ts = ts;}
    void convection_set(const ConvScheme & cs) {lim.set(cs);}

    Linear * solver;

    const Matter * fluid() const {return flu;}
    const Matter * solid() const {return sol;}

    const Domain * domain() const {return dom;}

  protected:
    TimeScheme diff_ts;
    TimeScheme conv_ts;
    Limiter    lim;

    const Domain * dom; 
    const Times  * time;
    Matter * flu;
    Matter * sol;

};

#endif
