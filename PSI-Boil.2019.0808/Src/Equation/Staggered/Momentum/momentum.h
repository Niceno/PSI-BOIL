#ifndef MOMENTUM_H
#define MOMENTUM_H

#include "../../../Parallel/mpi_macros.h"
#include <cmath>
#include <sstream>
#include <iostream>

#include "../staggered.h"
#include "../../../Matrix/matrix.h"
#include "../../../Parallel/Out/out.h"

/***************************************************************************//**
*  \brief Discretizes and solves momentum conservaion equation.
*
*  The eqation is discretized in integral form:                       
*  \f[
*        \int_{dV} \frac{\partial \rho {\bf u}}{\partial t} dV
*      + \int_{dS} \rho {\bf u} {\bf u} \, dS
*      = \int_{dS} \mu \nabla {\bf u} \, dS
       - \int_{dV} \nabla p \, dV
*      + {\bf F}
*      \; \; \; \;
*      [\frac{kg \, m}{s^2} = N]
*  \f] 
*  where \f$\rho \; [\frac{kg}{m^3}]\f$ is density, 
*  \f${\bf u} \; [\frac{m}{s}]\f$ is velocity, \f${t} \; [s]\f$ is time,
*  \f$\mu \; [\frac{kg}{m s}]\f$ is dynamic viscosity,
*  \f$p \; [\frac{kg}{m s^2}]\f$ is pressure  
*  \f${\bf F} \; [N]\f$ is a summ of external forces. 
*
*  Examples of external forces are gravity:
*  \f[
*      {\bf F}_g = \int_{dV} \rho {\bf g} \, dV
*      \; \; \; \; 
*      [\frac{kg}{m^3} \cdot \frac{m}{s^2} \cdot m^3 = \frac{kg \, m}{s^2}=N],
*  \f]
*  or pressure head:
*  \f[
*      {\bf F}_p = \int_{dV} \frac{\partial p}{\partial x} {\bf i} \; dV
*      \; \; \; \; 
*      [\frac{1}{m} \cdot \frac{kg}{m s^2} = \frac{kg \, m}{s^2} = N],
*  \f]
*  or surface tension force:
*  \f[
*      {\bf F}_s = \int_{dV} \sigma \kappa {\bf n} \; dV
*      \; \; \; \; 
*      [\frac{N}{m} \cdot \frac{1}{m^2} \cdot m^3 = N],
*  \f]
*  where \f$\sigma \; [\frac{N}{m}]\f$f is surface tension, property of the
*  fluids at interface and \f$ \kappa {\bf n} \; [\frac{1}{m^2}\f$ is the
*  surface curvature times unit surface normal. Details of the implementation
*  of surface tension force can be found in LevelSet.
*******************************************************************************/

////////////////
//            //
//  Momentum  //
//            //
////////////////
class Momentum : public Staggered {
	
  public:
    Momentum(const Vector & U, 
             const Vector & F, 
             Times & t, 
             Linear * sm,
             Matter * M);
    ~Momentum();

    /* interface to create_system */
    void discretize(const Scalar * mu_eddy = NULL) {
      boil::timer.start("momentum discretize");
      create_system(mu_eddy);
      boil::timer.stop("momentum discretize");
    }

    real cfl_max() const;
    void solve(const ResTol & toler, const ResRat & fact);
    void solve(const ResTol & toler) { solve(toler,ResRat(-1.)); };
    void solve(const ResRat & fact) { solve(ResTol(boil::atto),fact); };
    void solve_wo_outlet(const ResTol & toler, const ResRat & fact);
    void solve_wo_outlet(const ResTol & toler) { solve_wo_outlet(toler,ResRat(-1.)); };
    void solve_wo_outlet(const ResRat & fact) { solve_wo_outlet(ResTol(boil::atto),fact); };
    real bulk(const Comp & m, const real & coord) const;

    void get_eps(Scalar * src);
    void get_q(Scalar * src);
    void project(const Scalar & frc);
    void project(const Scalar & frc, Vector & veloc);
    void project_ghost(const Scalar & frc, const Scalar & c, const Scalar & k);
    void project_w_outlet(const Scalar & frc);
    void project_w_outlet(const Scalar & frc, Vector & veloc);
    void new_time_step(const Scalar * prs = NULL);
    void new_time_step(const Vector & v, const Scalar * prs = NULL);
    void grad(Scalar & p);
    void convection(const Scalar * prs = NULL) {convection(&cnew,prs);}
    real vol_phase_change(Scalar * psrc);
    
    void save(const char *, const int = -1);
    void load(const char *, const int = -1);

    void outlet();
    void pressure_outlet(const Scalar & frc);
    void pressure_outlet(const Scalar & frc, Vector & veloc);

    void vanishing_derivative_outlet(Vector & veloc);

    // trust velocities defined in solid (immersed boundary)
    void set_trust_vel_wall(const bool b){
           ib_trust_vel_wall=b; boil::oout<<"Momentum:trust_vel_wall= "<<b<<"\n";
        };
    bool get_trust_vel_wall(){return ib_trust_vel_wall;};

    void set_volf_ib(real x) {volf_ib=x;};

    void set_min_iteration(const int i){
      min_iter = MinIter(i);
      boil::oout<<"Momentum:min_iteration= "<<min_iter<<"\n";
    }
    int get_min_iteration() {return min_iter;}

    Matrix * A[3];

  private:
    
    void create_system(const Scalar * mu_eddy);

    void extrapolate_outlet_velocity(const real ubo, const real ratio);
    void convective_outlet(Vector & veloc, const real ubo);
    void scale_outlet_velocity(const real ubo, const real ratio);

    void convection(Vector * conv, const Scalar * prs = NULL);
    void diffusion();

    void advection_rho_c();
    void advection_rho_v();

    real bulk_i(const real & xp) const;
    real bulk_j(const real & yp) const;
    real bulk_k(const real & zp) const;

    real v_phase_change,volf_ib;
    bool ifull, jfull, kfull;
    bool ib_trust_vel_wall;
    MinIter min_iter;

};

#endif
