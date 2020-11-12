#ifndef AC_H
#define AC_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <cmath>
#include <array>
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
               const ResTol & toler, 
               const ResRat & factor, 
               const std::array<MaxIter,3> & mv,
               int * ncyc = NULL);

    bool vcycle(const ResTol & toler, const ResRat & factor, int * ncyc = NULL) {
      return cycle(Cycle::none(),Cycle::V(),toler,factor,mv_def,ncyc);
    }
    bool fcycle(const ResTol & toler, const ResRat & factor, int * ncyc = NULL) {
      return cycle(Cycle::none(),Cycle::F(),toler,factor,mv_def,ncyc);
    }
    bool wcycle(const ResTol & toler, const ResRat & factor, int * ncyc = NULL) {
      return cycle(Cycle::none(),Cycle::W(),toler,factor,mv_def,ncyc);
    }

    bool fullcycle(const Cycle & loop, const ResTol & toler, const ResRat & factor,
                   int * ncyc = NULL) {
      return cycle(Cycle::V(),loop,toler,factor,mv_def,ncyc);
    }

    bool vcycle(const ResRat & factor, int * ncyc = NULL) { 
      return cycle(Cycle::none(),Cycle::V(),ResTol(-1.0),factor,mv_def,ncyc);
    }
    bool fcycle(const ResRat & factor, int * ncyc = NULL) { 
      return cycle(Cycle::none(),Cycle::F(),ResTol(-1.0),factor,mv_def,ncyc);
    }
    bool wcycle(const ResRat & factor, int * ncyc = NULL) { 
      return cycle(Cycle::none(),Cycle::W(),ResTol(-1.0),factor,mv_def,ncyc);
    }

    bool fullcycle(const Cycle & loop, const ResRat & factor,
                   int * ncyc = NULL) { 
      return cycle(Cycle::V(),loop,ResTol(-1.0),factor,mv_def,ncyc);
    }

    bool vcycle(const ResTol & toler, int * ncyc = NULL) {
      return cycle(Cycle::none(),Cycle::V(),toler,ResRat(-1.0),mv_def,ncyc);
    }
    bool fcycle(const ResTol & toler, int * ncyc = NULL) {
      return cycle(Cycle::none(),Cycle::F(),toler,ResRat(-1.0),mv_def,ncyc);
    }
    bool wcycle(const ResTol & toler, int * ncyc = NULL) {
      return cycle(Cycle::none(),Cycle::W(),toler,ResRat(-1.0),mv_def,ncyc);
    }

    bool fullcycle(const Cycle & loop, const ResTol & toler,
                   int * ncyc = NULL) {
      return cycle(Cycle::V(),loop,toler,ResRat(-1.0),mv_def,ncyc);
    }

    //! Set and get the maximum number of cycles
    void max_cycles(const int mc) {max_cyc = mc;}
    int  max_cycles() const       {return max_cyc;}

    //! Set and get the minimum number of cycles
    void min_cycles(const int mc) {min_cyc = mc;}
    int  min_cycles() const       {return min_cyc;}

    //! Set and get max-iteration array
    void max_iter_array(const std::array<MaxIter,3> & mn) {mv_def = mn;}
    std::array<MaxIter,3> max_iter_array() const          {return mv_def;}

    //! Stop if diverging
    void stop_if_diverging(const bool & s) {stop_if_div = s;}
    bool stop_if_diverging() const         {return stop_if_div;}

    //! Use linf error
    void use_linf_error(const bool & s) {use_linf = s;}
    bool use_linf_error() const         {return use_linf;}

  private:
    //! Computes residual of the linear system. 
    real residual(Centered & h, real * linf = NULL) const;

    //! Restriction. Meaning from finer to coarser grid.
    void restriction(const Centered & h, Centered & H) const;
 
    //! Interpolation, Meaning from coarser to finer grid.
    void interpolation(const Centered & H, Centered & h) const;

    //! Coarsens immersed body flag.
    void coarsen_flag(const Centered & h, Centered & H) const; 

    //! Creates coarser discretized system.
    void coarsen_system(const Centered & h, Centered & H) const; 

    //! Individual components of a cycle
    bool init_cycles(const ResTol & toler, real & res_0, int * ncyc);
    int converged(const ResTol & toler, const ResRat & factor,
                  const int & cycle,
                  const real & res_0, const real & res0,
                  int * ncyc);
    void call_smoother(const int l, const MaxIter & mi,
                       const ResRat & res_rat, const ResTol & res_tol); 
    void call_solver(const int l, const MaxIter & mi,
                     const ResRat & res_rat, const ResTol & res_tol); 

    void v_cycle(const int l, const std::array<MaxIter,3> & mv);
    void f_cycle(const int l, const std::array<MaxIter,3> & mv);
    void w_cycle(const int l, const std::array<MaxIter,3> & mv);
    void full_cycle(const int l, const Cycle & cyc, const std::array<MaxIter,3> & mv);

    //! Pointers to coarser levels. 
    Centered * L[64];
    int        nlevels;

    //! Solver at coarsest level.
    Linear * solver;

    //! Parameters for steering of a cycle
    int  max_cyc;
    int  min_cyc;
    bool stop_if_div, use_linf;
    std::array<MaxIter,3> mv_def;
};

#endif
