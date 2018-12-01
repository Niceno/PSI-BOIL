#ifndef ENTHALPYTIF_H
#define ENTHALPYTIF_H

#include "../../../Parallel/mpi_macros.h"
#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Solver/Gauss/gauss.h"
#include "../../../Timer/timer.h"

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
        \param f   - extarnal source array (\f$\dot{q}\f$),
        \param u   - convection velocity (\f${\bf u}\f$),
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - Krylov subspace solver. It acts as a solver, or as a
                     smoother for AC multirid.
        \param flu - Holds all fluid properties (\f$\rho, C_p, \lambda\f$),
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
        \param tsat - (reference) saturation temperature.
        \param latent - latent heat of evaporation.
        \param mresis - interfacial mass transfer resistance.
        \param fs - free surface position from VOF.
    */

    /* Most general constructor */
    EnthalpyTIF(const Scalar & phi, 
               const Scalar & f,
               const Scalar & clr,
               const Vector & u,
               Times & t,
               Krylov * sm,
               Matter * flu,
               const real tsat,
               Matter * sol = NULL,
               const Vector * fs = NULL,
               const real latent = 1.0,
               const real mresis = 0.0,
               const Scalar * mdot = NULL,
               const Scalar * pres = NULL);

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
                Matter * sol = NULL) :
      EnthalpyTIF(phi,f,clr,u,t,sm,flu,tsat,sol) { }
#endif
    ~EnthalpyTIF();

    void new_time_step(const Scalar * diff_eddy = NULL);
    void solve(const ResRat & fact, const char * name = NULL);
    void solve_sor(const int & it, const real & r, const char * name = NULL);
    void deltat(Scalar & deltaT, const Scalar & heaviside, 
                const ResRat & fact, const char * name = NULL, const real flag = blendfactor);
    void update_tifold() {
      for_ijk(i,j,k) 
        tifold[i][j][k] = tif[i][j][k];
      //boil::oout<<"EnthalpyFD::tint_field()  initialize tif"<<"\n";
      store_tif = true;
      tifold.exchange();
    }

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

    void tint_field(const Scalar & heaviside, const real factor = blendfactor, const bool iter = false);
    Scalar tif, tifold;

    real get_blendfactor(){return blendfactor;}
    void set_blendfactor(real b){
      blendfactor=b;
      boil::oout<<"EnthalpyTIF:blendfactor= "<<blendfactor<<"\n";
    }

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
    real tsat,rhol,rhov,cpl,cpv,lambdal,lambdav,clrsurf,epsl;
    bool store_clrold;
    Scalar clrold;
    Scalar ftif,ftifold,fdelta;
    ScalarInt iflag;
    real turbP; /* turbulent Prandtl number */
    bool laminar;

    const Vector * fs;
    Vector fsold;

    real latent, mresis;
    static real blendfactor;
    bool store_tif;
    const Scalar * mflx;
    const Scalar * pres;
    void Pressure_effect(const Scalar & heaviside);
    void Mass_src_effect(const Scalar & heaviside);
    void Extend_tint    (const Scalar & heaviside);
    void update_ftif(const Scalar * diff_eddy = NULL);

    bool Vicinity(const int i, const int j, const int k,
                  const Scalar & heaviside);
    //bool Vicinity(const int i, const int j, const int k);

    bool Interface(const int i, const int j, const int k,
                   const Scalar & heaviside);
    bool Interface(const real heavi);

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

