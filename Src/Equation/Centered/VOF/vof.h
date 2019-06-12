#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_realistic.h"
#include "../../Heaviside/heaviside.h"
#include "../../Topology/topology.h"

///////////
//       //
//  VOF  //
//       //
///////////
class VOF : public Centered {
  public:
    VOF(const Scalar & phi,
        const Scalar & f,
        const Scalar & kappa,
        const Vector & u, 
        Times & t,
        Krylov * S,
        Vector * bndclr = NULL,
        Matter * flu = NULL);
    ~VOF();

    void new_time_step(){};
    void advance();
    void advance(Scalar & sca);
    void curvature();
    void ancillary(); /* calcs ancillary params such as adens w/o advance */
    void tension(Vector * vec, const Matter matt);
    void totalvol();
    void front_minmax();
    void init(){};
    void smooth(const Scalar & sca, Scalar & scb, const int itnum);

    // getter for front_minmax
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};

    // getter and setter for wall value tolerance
    real get_tol_wall() { return tol_wall; }
    void set_tol_wall(real tolnew) {
      tol_wall = tolnew;
      boil::oout<<"VOF: New wall value tolerance: "<<tol_wall<<boil::endl;
      return;
    }

    // getter and setter for curv_method
    int get_curv_method() {return curv_method;}
    void set_curv_method(int i) {
      curv_method=i;
      if(i==0){
        boil::oout<<"VOF: height function is used for curvature calculation.\n";
      } else if(i==1){
        boil::oout<<"VOF: smoothed VOF is used for curvature calculation.\n";
      } else {
        boil::oout<<"method should be 0 or 1.\n";
        boil::oout<<"0 for height function.\n";
        boil::oout<<"1 for smoothed VOF.\n";
        exit(0);
      }
    }

    /* setter for cangle */
    void set_cangle(const real r) {
      cangle=r/180.0*acos(-1.0);
      boil::oout<<"set_cangle: cangle= "<<r<<"\n";
      boil::oout<<"#############################################\n";
      boil::oout<<"# WARNING!!! cangle is not implemented yet. #\n";
      boil::oout<<"#############################################\n";
    }
    /* getter for cangle */
    real get_cangle() { return(cangle/acos(-1.0)*180.0);};

    Vector fs;
    Vector * bndclr;
    Topology  topo;

    Scalar nalpha;
    Scalar nx,ny,nz;/* normal to interface */
    Scalar adens;
    Scalar mx,my,mz;/* normal to interface, in real space */

  protected:
    void advance_x(Scalar & sca);
    void advance_y(Scalar & sca);
    void advance_z(Scalar & sca);
    void bdcurv(const Scalar & g, const real & v);
    void cal_fs3();
    void cal_fs_interp();
    void fs_bnd();
    void curv_HF();
    void curv_smooth();
    void extract_alpha();
    void true_norm_vect();
    void standardized_norm_vect();
    void insert_bc(const Scalar & g);
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_norm_cc(const Scalar & g);
    void insert_bc_norm();
    void norm_cc(const Scalar & g);
    void norm_young(const Scalar & g);
    void normalize(real & r1, real & r2, real & r3);
    void wall_adhesion_norm(real & nx, real & ny, real & nz,
                            const real nwx, const real nwy, const real nwz);
    void update_at_walls();
    real kappa_ave(const real r1, const real r2);
    real kappa_ave(const real r1, const real r2, const int i1, const int i2);

    real calc_v(real r1, real r2, real r3, real r4);
    real calc_alpha(const real r1, const real r2, const real r3, const real r4);
    real calc_flux(const real g, real c, const real nx, const real ny, const real nz);

    void selectMax(const real r1, const real r2, const real r3,
                   const real r4, const real r5, const real r6,
                   const real r7, const real r8, const real r9,
                   const int i1,  const int i2,  const int i3);
    void set_iflag();
    void insert_bc_flag(ScalarInt & g, const bool b);

    void cal_adens();
    void cal_bndclr();

    void set_adens(const Scalar & newadens) {
      for_aijk(i,j,k)
        adens[i][j][k] = newadens[i][j][k];
    }

    real extrapolate_v(const int i, const int j, const int k,
                       const int ofx, const int ofy, const int ofz,
                       const real xp, const real yp, const real zp);

    void norm_cc_imin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_imax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmax(const Scalar &g, const int i,const int j, const int k);
    void vf_limiter();   

    real alpha_val(const int i, const int j, const int k);
    real fs_val(const Comp m, const int i, const int j, const int k);
    real frontPosition(const int i, const int j, const int k, const Comp m);

    Scalar clr;     /* color function */
    Scalar kappa;        /* curvature */
    Scalar stmp;
    ScalarInt iflag,iflagx,iflagy,iflagz;

    real rhol, rhov; /* densities for velocity correction */
    const Matter * mixt() const {return mixture;}
    Matter * mixture;

    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real theta;
    real epsnorm;
    real phisurf;
    real tol_wall, tol_flux, tol_ext, flux_cfl;
    real ww, dxmin;
    bool iminp, imaxp, jminp, jmaxp, kminp, kmaxp; // true = periodic
    bool iminw, imaxw, jminw, jmaxw, kminw, kmaxw; // true = wall
    bool iminc, imaxc, jminc, jmaxc, kminc, kmaxc; // true = cut-stencil

    Heaviside heavi;

    int nlayer;
    int curv_method;
    real cangle;
};	
#endif

