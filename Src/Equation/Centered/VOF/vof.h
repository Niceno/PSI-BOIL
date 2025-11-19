#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_realistic.h"
#include "../../Heaviside/heaviside.h"
#include "../../Topology/topology.h"

#define IB

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
        Vector * bndclr = NULL);
    ~VOF();

    void new_time_step(){};
    void advance(const bool anci = true);
    void advance(Scalar & sca, const bool anci = true);
    void curvature();
    void ancillary(); /* calcs ancillary params such as adens w/o advance */
    void ancillary(Scalar & scp);
    void tension(Vector * vec, const Matter matt);
    void tension(Vector * vec, const Matter matt, const Scalar & scp);
    void totalvol();
    void front_minmax();
    void front_minmax( Range<real> xr
                     , Range<real> yr
                     , Range<real> zr );

    void init(){ ancillary(); };

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
    }
    /* getter for cangle */
    real get_cangle() { return(cangle/acos(-1.0)*180.0);};

    /* setter for limit_color */
    void set_limit_color(const bool b) {
      limit_color=b;
      boil::oout<<"set_limit_color= "<<b<<"\n";
    }
    /* getter for limit_color */
    bool get_limit_color() { return(limit_color);};

    /* setter for use_subgrid */
    void set_use_subgrid(const bool b) {
      use_subgrid=b;
      boil::oout<<"set_use_subgrid= "<<b<<"\n";
    }
    /* getter for use_subgrid */
    bool get_use_subgrid() { return(use_subgrid);};

    /* setter for use_HF_wall */
    void set_use_HF_wall(const bool b) {
      use_HF_wall=b;
      if (use_HF_wall) {
        boil::oout<<"##################### Warning ########################\n";
        boil::oout<<"## Height function is used for wall adhesion force. ##\n";
        boil::oout<<"## Equilibrium contact angle must be 0 or 180       ##\n";
        boil::oout<<"######################################################\n";
      }
    }
    /* getter for use_subgrid */
    bool get_use_HF_wall() { return(use_HF_wall);};

    /* min and max of color function in fluid domain */
    real minval() {return minclr;}
    real maxval() {return maxclr;}
    void color_minmax();
    void set_minval(real r) {minclr=r;}
    void set_maxval(real r) {maxclr=r;}

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
    void bdcurv();
    void cal_fs3(const Scalar & scp);
    void cal_fs_interp(const Scalar & scp);
    void curv_HF();
    void curv_smooth();
    void extract_alpha(Scalar & scp);
    void extract_alpha_near_bnd(const Scalar & scp);
    void fs_bnd(const Scalar & scp);
    void fs_bnd_nosubgrid(const Scalar & scp);
    void standardized_norm_vect();
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void ib_norm(const Scalar & g);
    void ib_norm_cal(const int cc, const int i, const int j, const int k);
    void insert_bc(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_gradphi(const Scalar & g);
    void insert_bc_norm_cc(const Scalar & g);
    void insert_bc_norm();
    void nib(const real & n1, const real & n2, const real & n3
           , const real & n4, const real & n5, const real & n6
           , real r[]);

    void norm_cc(const Scalar & g);
    void norm_cc(real & nx_val, real & ny_val, real & nz_val, Comp & mcomp,
                 const int i, const int j, const int z, const Scalar & sca);

    void norm_young(const Scalar & g, const bool b);
    void norm_young(real & nx_val, real & ny_val, real & nz_val,
                    const int i, const int j, const int z, const Scalar & sca,
                    const bool b);

    void norm_mixed(const Scalar & g);
    void norm_mixed(real & nx_val, real & ny_val, real & nz_val,
                    const int i, const int j, const int k,
                    const Scalar & sca);

    void bdnorm(Scalar & scp);
    void extend_norm(const Scalar & g);

    void normalize(real & r1, real & r2, real & r3);
    void normalize_l1(real & nx_l1, real & ny_l1, real & nz_l1,
                      const real nx, const real ny, const real nz);

    void nwall(const Scalar & g
             , const real & r1, const real & r2, const real & r3
             , const int & i1, const int & i2, const int &i3
             , real r[] );
    void smooth(const Scalar & sca, Scalar & scb, const int itnum);
    void true_norm_vect();
    void update_at_walls(Scalar & scp, const real * cang = nullptr);
    void wall_norm(const Scalar & sca);
    real kappa_ave(const real r1, const real r2);
    real kappa_ave(const real r1, const real r2, const int i1, const int i2);

    real calc_v(real r1, real r2, real r3, real r4);
    real calc_alpha(const real r1, const real r2, const real r3, const real r4);
    real calc_flux(const real g, real c, const real nx, const real ny, const real nz);

    void select_norm_cc(real & nx_val, real & ny_val, real & nz_val,
                        real & NxX, real & NyX, real & NzX,
                        real & NxY, real & NyY, real & NzY,
                        real & NxZ, real & NyZ, real & NzZ,
                        Comp & mcomp);
    void select_norm_myc(real & nx_val, real & ny_val, real & nz_val,
                         const real & nx_cc, const real & ny_cc, const real & nz_cc,
                         const real & nx_young, const real & ny_young, const real & nz_young,
                         const Comp & mcomp);


    void set_iflag();
    void insert_bc_flag(ScalarInt & g, const bool b);

    void cal_adens();
    void cal_bndclr(const Scalar & scp);

    void set_adens(const Scalar & newadens) {
      for_aijk(i,j,k)
        adens[i][j][k] = newadens[i][j][k];
    }

    real extrapolate_v(const int i, const int j, const int k,
                       const int ofx, const int ofy, const int ofz,
                       const real xp, const real yp, const real zp,
                       const Scalar & scp, const real * cang = nullptr);

    void norm_cc_imin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_imax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmax(const Scalar &g, const int i,const int j, const int k);
    void vf_limiter();   

    /* elvira functions */
    void norm_elvira(Scalar & sca);
    void norm_elvira(int i,int j,int k,
                     real valcc,real valmc,real valpc,real valcm,
                     real valcp,real valmm,real valpm,real valmp,real valpp);
    void normalize_elvira(const real m, const real sig, 
                          real & nnx, real & nny, real & nnz);
    real elvira_l2(const real alp,const real nnx,const real nny,const real nnz,
                   const real valcc,const real valmc,const real valpc,
                   const real valcm,const real valcp,const real valmm,
                   const real valpm,const real valmp,const real valpp);
    real ext_v(const real xp, const real yp, const real zp, 
               const real vv1, const real vv2, const real vv3,
               const real vn1, const real vn2, const real vn3,
               const real denom, const real alp);

    real alpha_val(const real c, const real nnx, const real nny, const real nnz);
    real fs_val(const Comp m, const int i, const int j, const int k);
    real frontPosition(const int i, const int j, const int k, const Comp m);

    Scalar kappa;        /* curvature */
    Scalar stmp,stmp2;
    ScalarInt iflag,iflagx;

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
    bool ifull, jfull, kfull; // true = not a dummy direction
    bool limit_color, use_subgrid, use_HF_wall;
    real minclr, maxclr;

    Heaviside heavi;

    int nlayer;
    int curv_method;
    real cangle;
};	
#endif

