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

//////////
//      //
//  AC  //
//      //
//////////
class AC {
	
  public:
    //! Basic constructor
    AC(Centered * cen); 

    //! The solution algorythm. V-cycle. 
    bool vcycle(const ResRat & factor, int * ncyc = NULL);

    //! Another way to call vcycle (temporary) 
    bool vcycle(int * ncyc = NULL) {return vcycle(ResRat(targ_res_rat), ncyc);}

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

    //! Pointers to coarser levels. 
    Centered * L[64];
    int        nlevels;

    //! Parameters for steering of v-cycle
    int  max_cyc;
    int  min_cyc;
    bool stop_if_div;
    real targ_res_val;
    real targ_res_rat;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: additive.h,v 1.20 2014/08/06 07:43:11 sato Exp $'/
+-----------------------------------------------------------------------------*/
