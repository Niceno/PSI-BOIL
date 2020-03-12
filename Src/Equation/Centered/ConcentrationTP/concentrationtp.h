#ifndef CONCENTRATIONTP_H
#define CONCENTRATIONTP_H

#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../Heaviside/MarchingSquares/MSaxisym/ms_axisym.h"
#include "../../Heaviside/MarchingCubes/marching_cubes.h"

/***************************************************************************//**
*  \brief Discretizes and solves species conservaiton equation in gas phase.
*
*  The eqation is discretized in integral form:
*  \f[
*        \int_V \frac{\partial \rho \varepsilon}{\partial t} dV
*      + \int_S \rho {\bf u_{vol}} \varepsilon \, dS
*      = \int_S \gamma \nabla \varepsilon \, dS
*      + \dot{M}
*      \; \; \; \;
*      [\frac{kg}{s}]
*  \f]
*  where \f$\varepsilon \; [1]\f$ is species volume fraction from 0 to 1,
*  \f$\rho \; [\frac{kg}{m^3}]\f$ is speciesmass density, 
*  \f${t} \; [s]\f$ is time,
*  \f${\bf u_{vol}} \; [\frac{m}{s}]\f$ is volumetric velocity,
*  \f$\gamma \; [\frac{kg}{ms}]\f$ is species diffusivity and
*  \f$\dot{M} \; [\frac{kg}{s}]\f$ is (external) mass source rate.
*******************************************************************************/

/////////////////////
//                 //
//  Concentration  //
//                 //
/////////////////////
class ConcentrationTP : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - species concentration array (\f$\epsilon\f$),
        \param f   - external source array (\f$\dot{m}\f$),
        \param u   - volume averaged velocity (\f${\bf u}\f$),
        \param color  - volume fraction array (\f$\phi\f$),
        \param flxclr - volume fraction flow array (\f$\phi\f$),
        \param heavi  - phase indicator object
        \param topo   - topology object
        \param t      - simulation (physical) time (\f${t}\f$),
        \param sm   - Linear solver. It acts as a solver, or as a
                      smoother for AC multigrid.
        \param flu  - Holds material properties of diffusing species.
    */
    ConcentrationTP(const Scalar & phi,
                    const Scalar & f,
                    const Vector & u,
                    const Scalar & color, 
                    const Vector & flxclr,
                    Heaviside * heavi,
                    Topology * topo,
                    Times & t,
                    Linear * sm,
                    Matter * flu); 
    ~ConcentrationTP();
	  
    //! Interface call to parent's discretization.
    void discretize(const Scalar * eddy = NULL);

    //! Start a new time step.
    void new_time_step(const Scalar * diff_eddy = NULL);
    void solve(const ResRat & fact, const char * name = NULL);

    inline real get_turbS() const { return turbS; }
    void set_turbS(const real a) {
      turbS=a;
      boil::oout<<"ConcentrationTP::turbS= "<<turbS<<"\n";
    }

    void convection();
    void extrapolate();

    //! Velocity calculations.
    void compute_umass(Vector & umass, const Vector & uvol,
                       const Matter * gas, const Scalar * diff_eddy = NULL);
    void compute_udiff_div(Scalar & p, const Property * mu_fluid,
                           const Vector & udiff);

    void init_color(const Scalar &clro) {
      for_aijk(i,j,k) clrold[i][j][k] = clro[i][j][k];
      store_clrold = true;
      boil::oout << "ConcentrationTP::initialize clrold\n";
    }

    void store() {
      /* store clrold */
      for_aijk(i,j,k){
        clrold[i][j][k] = (*clr)[i][j][k];
      }
    }

  protected:
    void create_system_innertial();
    void create_system_diffusive(const Scalar * diff_eddy = NULL);

    void convection(Scalar * sca);

    /* color function on faces */
    real surface_color(const Comp & m, const int i, const int j, const int k);
    real surface_color(const Sign sig, const Comp & m,
                       const int i, const int j, const int k);

    /* extrapolation */
    void extrapolation_flag();

    const Scalar * clr;
    Heaviside * heavi;
    Topology * topo;
    Scalar clrold;
    ScalarInt eflag;
    const Vector * colorflow;
    const Property * rho_dif;
    const Property * dcoef;

    bool store_clrold;
    real turbS; /* turbulent schmidt number */
    bool laminar;
};	
#endif