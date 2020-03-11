#ifndef AC_H
#define AC_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <cmath>
#include "../../Equation/Centered/centered.h"
#include "../../Timer/timer.h"

/***************************************************************************//**
*  \brief Implementation of the Additive Correction (AC) multigrid solver.
*
*  Implements the AC multigrid solver, described by Raithby (...).
*  It is defined for Centered variables, such as Pressure, LevelSet, 
*  Temperature, Concentration ...
*******************************************************************************/

#include "additive_ravioli.h"

//////////
//      //
//  AC  //
//      //
//////////
class AC {
	
  public:
    //! Basic constructor
    AC(Centered * cen, Linear * sol = NULL); 

    //! The solution algorithm + shorthand calls added for simplicity. 
    bool cycle(const Cycle & init, const Cycle & loop, 
               const ResRat & factor, int * ncyc = NULL);

    bool vcycle(const ResRat & factor, int * ncyc = NULL) { 
      return cycle(Cycle::none(),Cycle::V(),factor,ncyc);
    }
    bool fcycle(const ResRat & factor, int * ncyc = NULL) { 
      return cycle(Cycle::none(),Cycle::F(),factor,ncyc);
    }
    bool wcycle(const ResRat & factor, int * ncyc = NULL) { 
      return cycle(Cycle::none(),Cycle::W(),factor,ncyc);
    }

    bool fullcycle(const Cycle & loop, const ResRat & factor,
                   int * ncyc = NULL) { 
      return cycle(Cycle::V(),loop,factor,ncyc);
    }

    //! Another way to call cycles (temporary) 
    bool cycle(const Cycle & init, const Cycle & loop,
               int * ncyc = NULL) {
      return cycle(init,loop,ResRat(targ_res_rat), ncyc);
    }

    bool vcycle(int * ncyc = NULL) {return vcycle(ResRat(targ_res_rat), ncyc);}
    bool fcycle(int * ncyc = NULL) {return fcycle(ResRat(targ_res_rat), ncyc);}
    bool wcycle(int * ncyc = NULL) {return wcycle(ResRat(targ_res_rat), ncyc);}

    bool fullcycle(const Cycle & loop, int * ncyc = NULL) {
      return fullcycle(loop,ResRat(targ_res_rat), ncyc);
    }

    //! Set and get the maximum number of cycles
    void max_cycles(const int mc) {max_cyc = mc;}
    int  max_cycles() const       {return max_cyc;}

    //! Set and get the minimum number of cycles
    void min_cycles(const int mc) {min_cyc = mc;}
    int  min_cycles() const       {return min_cyc;}

    //! Set and get target residual value
    void target_residual_value(const real & r) {targ_res_val = r;}
    real target_residual_value() const         {return targ_res_val;}

    //! Set and get target residual ratio
    void target_residual_ratio(const real & r) {targ_res_rat = r;}
    real target_residual_ratio() const         {return targ_res_rat;}

    //! Stop id diverging
    void stop_if_diverging(const bool & s) {stop_if_div = s;}
    bool stop_if_diverging() const         {return stop_if_div;}

  private:
    //! Computes residual of the linear system. 
    real residual(Centered & h) const;

    //! Restriction. Meaning from finer to coarser grid.
    void restriction(const Centered & h, Centered & H) const;
 
    //! Interpolation, Meaning from coarser to finer grid.
    void interpolation(const Centered & H, Centered & h) const;

    //! Creates coarser discretized system
    void coarsen_system(const Centered & h, Centered & H) const; 

    //! Individual components of a cycle
    bool init_cycles(real & res_0, int * ncyc);
    int converged(const ResRat & factor, const int & cycle,
                  const real & res_0, const real & res0,
                  int * ncyc);
    void call_smoother(const int l, const MaxIter & mi,
                       const ResRat & res_rat, const ResTol & res_tol); 
    void call_solver(const int l, const MaxIter & mi,
                     const ResRat & res_rat, const ResTol & res_tol); 

    void v_cycle(const int l);
    void f_cycle(const int l, const bool upstream);
    void w_cycle(const int l, const bool upstream);
    void full_cycle(const int l, const Cycle & cyc);

    //! Pointers to coarser levels. 
    Centered * L[64];
    int        nlevels;

    //! Solver at coarsest level.
    Linear * solver;

    //! Parameters for steering of a cycle
    int  max_cyc;
    int  min_cyc;
    bool stop_if_div;
    real targ_res_val;
    real targ_res_rat;
};

#endif
