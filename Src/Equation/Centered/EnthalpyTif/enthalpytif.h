#ifndef ENTHALPYTIF_H
#define ENTHALPYTIF_H

#include "../../../Parallel/mpi_macros.h"
#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Solver/Gauss/gauss.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../Tifmodel/tif.h"

/***************************************************************************//**
*  \brief Discretizes and solves enthalpy conservaion equation.
*
*  The eqation is discretized in integral form:                       
*  \f[
*        \int_V \frac{\partial \rho C_p T}{\partial t} dV
*      + \int_S \rho C_p {\bf u} T \, dS
*      = \int_S \lambda \nabla T \, dS
*      + \dot{Q}
*      \; \; \; \;
*      [\frac{J}{s} = W]
*  \f] 
*  where \f$T \; [K]\f$ is temperature, \f$\rho \; [\frac{kg}{m^3}]\f$ is 
*  density, \f$C_p \; [\frac{J}{kgK}]\f$ is thermal capacity, \f${t} \; [s]\f$ 
*  is time, \f${\bf u} \; [\frac{m}{s}]\f$ is convective velocity,
*  \f$\lambda \; [\frac{W}{mK}]\f$ is thermal conductivity and 
*  \f$\dot{Q} \; [\frac{J}{s}]\f$ is (external) heat source rate. 
*******************************************************************************/

////////////////
//            //
//  Enthalpy  //
//            //
////////////////
class EnthalpyTIF : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - temperature (\f$T\f$),
        \param f   - external source array (\f$\dot{q}\f$),
        \param u   - convection velocity (\f${\bf u}\f$),
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - Krylov subspace solver. It acts as a solver, or as a
                     smoother for AC multirid.
        \param flu - Holds all fluid properties (\f$\rho, C_p, \lambda\f$),
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
        \param tifmodel - interfacial temperature model.
        \param latent - latent heat of evaporation.
        \param mresis - interfacial mass transfer resistance.
        \param fs - free surface position from VOF.
        \param adens - interfacial area density.
    */

    /* Most general constructor */
    EnthalpyTIF(const Scalar & phi, 
                const Scalar & f,
                const Scalar & clr,
                const Vector & u,
                Times & t,
                Krylov * sm,
                Matter * flu,
                TIF & tifmodel,
                Matter * sol = NULL,
                const Vector * fs = NULL,
                const Scalar * adens = NULL);

    /* would it be possible to construct an ETIF object using real tsat 
     * for bckwrd compatibility? */

#if 0  /* sadly, constructors calling different constructors requires C++11 */
    /* Original CIPCSL2-EnthalpyFD constructor */
    EnthalpyTIF(const Scalar & phi,
                const Scalar & f,
                const Scalar & clr,
                const Vector & u,
                Times & t,
                Krylov * sm,
                Matter * flu,
                const real tsat,
                Matter * sol = NULL);
#endif
    ~EnthalpyTIF();

    void new_time_step(const Scalar * diff_eddy = NULL);
    void solve(const ResRat & fact, const char * name = NULL);
    void solve_sor(const int & it, const real & r, const char * name = NULL);

#if 0
    //! Direct solver introduced just for checking it.
    void direct() {
      for_ijk(i,j,k)
        fnew[i][j][k] = fold[i][j][k]
                      + cnew[i][j][k] * conv_ts.N()
                      + fbnd[i][j][k]
                      + fext[i][j][k];

      boil::timer.start("enthalpy solver");
      Gauss gs;
      gs.solve(A, phi, fnew);
      boil::timer.stop("enthalpy solver");
    }
#endif

    real hflux_wall(const Scalar & s, const Dir d
                  , const Scalar * diff_eddy = NULL);
    real hflux_wall_ib(const Scalar * diff_eddy = NULL);

    //! Interface call to parent's discretization.
    void discretize(const Scalar * diff_eddy = NULL) {
      boil::timer.start("enthalpy discretize");
      create_system(diff_eddy);
      if(diff_eddy != NULL)laminar=false;
      boil::timer.stop("enthalpy discretize");
    }

    real get_turbP(){return turbP;}
    void set_turbP(real a){
      turbP=a;
      boil::oout<<"EnthalpyTIF:turbP= "<<turbP<<"\n";
    }

    void convection();

    /* to be removed */
    void update_ftif(const Scalar * diff_eddy = NULL);

  protected:
    void create_system(const Scalar * diff_eddy = NULL);
    void create_system_innertial();
    void create_system_diffusive(const Scalar * diff_eddy = NULL);
    void create_system_bnd();
    real update_rhs();
    void convection(Scalar * sca);
    void diffusion_fd(const Scalar * diff_eddy = NULL);
    real dVFD(const int i, const int j, const int k){
      real vol = 0.5 * (phi.dxw(i)+phi.dxe(i))
               * 0.5 * (phi.dys(j)+phi.dyn(j))
               * 0.5 * (phi.dzb(k)+phi.dzt(k));
      return vol;
    };
    void setflag();
    void diff_matrix(real & am, real & ac, real & ap
                , real & tm, real & tc, real & tp
                , real & aflagm, real & aflagp
                , const real vol, const real area
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const real clm, const real clc, const real clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , real pm, real pc, real pp
                , const real edm, const real edc, const real edp
                , const int i, const int j, const int k, const Comp m);

    const Scalar * clr;
    real rhol,rhov,cpl,cpv,lambdal,lambdav,clrsurf,epsl;
    bool store_clrold;
    Scalar clrold;
    Scalar ftif,ftifold; /* tbd */
    ScalarInt iflag;
    real turbP; /* turbulent Prandtl number */
    bool laminar;

    const Scalar * adens;
    Scalar adensold;

    const Vector * fs;
    Vector fsold;

    TIF & tifmodel;  

    bool Interface(const int dir, const Comp m,
                   const int i, const int j, const int k);

    real Tint(const int dir, const Comp &mcomp, const real frac,
              const int i, const int j, const int k);
    real Tint_old(const int dir, const Comp &mcomp, const real frac,
              const int i, const int j, const int k);
    real Tint(const int i, const int j, const int k);
    real Tint_old(const int i, const int j, const int k);

    real distance_x(const int i, const int j, const int k,
                    const int dir, real & tint,
                    const bool old = false);
    real distance_y(const int i, const int j, const int k,
                    const int dir, real & tint,
                    const bool old = false);
    real distance_z(const int i, const int j, const int k,
                    const int dir, real & tint,
                    const bool old = false);
    bool distance1D_x(const int i, const int j, const int k,
                      const int dir, real & tint, real & dist);
    bool distance1D_y(const int i, const int j, const int k,
                      const int dir, real & tint, real & dist);
    bool distance1D_z(const int i, const int j, const int k,
                      const int dir, real & tint, real & dist);
    bool distance1D_xold(const int i, const int j, const int k,
                         const int dir, real & tint, real & dist);
    bool distance1D_yold(const int i, const int j, const int k,
                         const int dir, real & tint, real & dist);
    bool distance1D_zold(const int i, const int j, const int k,
                         const int dir, real & tint, real & dist);
};	
#endif

