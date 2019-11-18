#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_realistic.h"
#include "../../Heaviside/heaviside.h"
#include "../../Topology/topology.h"

#define IB

using arr3D = std::vector< std::vector< std::vector<real> > >;

#include "vof_ravioli.h"

///////////
///      //
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
    void forward(Scalar & scp);
    virtual void advance(const bool anci = true);
    virtual void advance(Scalar & sca, const bool anci = true);
    void curvature();
    virtual void ancillary(); /* calcs ancillary params such as adens w/o advance */
    void tension(Vector * vec, const Matter matt);
    void tension(Vector * vec, const Matter matt, const Scalar & scp);
    void totalvol();
    void front_minmax();
    void front_minmax( Range<real> xr
                     , Range<real> yr
                     , Range<real> zr );

    void init(){ ancillary(); };

    // getter for front_minmax
    inline real get_xminft() { return(xminft);};
    inline real get_xmaxft() { return(xmaxft);};
    inline real get_yminft() { return(yminft);};
    inline real get_ymaxft() { return(ymaxft);};
    inline real get_zminft() { return(zminft);};
    inline real get_zmaxft() { return(zmaxft);};

    // getter and setter for wall value tolerance
    inline real get_tol_wall() { return tol_wall; }
    inline void set_tol_wall(real tolnew) {
      tol_wall = tolnew;
      boil::oout<<"VOF: New wall value tolerance: "<<tol_wall<<boil::endl;
      return;
    }

    // getter and setter for bulk_curv_method
    inline CurvMethod get_curv_method() {return bulk_curv_method;}
    void set_curv_method(const CurvMethod cm) {
      bulk_curv_method=cm;
      if(cm==CurvMethod::HF()){
        boil::oout<<"VOF: height function is used for curvature calculation.\n";
      } else if(cm==CurvMethod::DivNorm()){
        boil::oout<<"VOF: smoothed VOF is used for curvature calculation.\n";
      } else {
        boil::oout<<"Curv method should be HF or DivNorm.\n";
        boil::oout<<"Exiting."<<boil::endl;
        exit(0);
      }
    }

    /* setter for cangle */
    inline void set_cangle(const real r) {
      cangle=r/180.0*boil::pi;
      boil::oout<<"set_cangle: cangle= "<<r<<"\n";
    }
    /* getter for cangle */
    inline real get_cangle() { return(cangle/boil::pi*180.0);};

    /* setter for limit_color */
    inline void set_limit_color(const bool b) {
      limit_color=b;
      boil::oout<<"set_limit_color= "<<b<<"\n";
    }
    /* getter for limit_color */
    inline bool get_limit_color() { return(limit_color);};

    /* setter for use_subgrid */
    inline void set_use_subgrid(const bool b) {
      use_subgrid=b;
      boil::oout<<"set_use_subgrid= "<<b<<"\n";
    }
    /* getter for use_subgrid */
    inline bool get_use_subgrid() { return(use_subgrid);};

    /* setter for use_interp */
    inline void set_use_interp(const bool b) {
      use_interp=b;
      boil::oout<<"set_use_interp= "<<b<<"\n";
    }
    /* getter for use_interp */
    inline bool get_use_interp() { return(use_interp);};

    /* min and max of color function in fluid domain */
    inline real minval() {return minclr;}
    inline real maxval() {return maxclr;}
    void color_minmax();
    inline void set_minval(real r) {minclr=r;}
    inline void set_maxval(real r) {maxclr=r;}

    /* setter for normal vector method */
    void set_normal_vector_method_advance(const NormMethod nm) {
      norm_method_advance = nm;
      boil::oout<<"Normal vector method for advance: "<<nm<<boil::endl;
      if(nm==NormMethod::ElviraYZ()) {
        mcomp_for_elvira = Comp::i();
      } else if(nm==NormMethod::ElviraXZ()) {
        mcomp_for_elvira = Comp::j();
      } else if(nm==NormMethod::ElviraXY()) {
        mcomp_for_elvira = Comp::k();
      } else {
        mcomp_for_elvira = Comp::undefined();
      }
    }
    void set_normal_vector_method_curvature(const NormMethod nm) {
      norm_method_curvature = nm;
      boil::oout<<"Normal vector method for curvature: "<<nm<<boil::endl;
      if(nm==NormMethod::ElviraYZ()) {
        mcomp_for_elvira = Comp::i();
      } else if(nm==NormMethod::ElviraXZ()) {
        mcomp_for_elvira = Comp::j();
      } else if(nm==NormMethod::ElviraXY()) {
        mcomp_for_elvira = Comp::k();
      } else {
        mcomp_for_elvira = Comp::undefined();
      }
    }
    void set_normal_vector_method_all(const NormMethod nm) {
      norm_method_advance = nm;
      norm_method_curvature = nm;
      boil::oout<<"Normal vector method: "<<nm<<boil::endl;
      if(nm==NormMethod::ElviraYZ()) {
        mcomp_for_elvira = Comp::i();
      } else if(nm==NormMethod::ElviraXZ()) {
        mcomp_for_elvira = Comp::j();
      } else if(nm==NormMethod::ElviraXY()) {
        mcomp_for_elvira = Comp::k();
      } else {
        mcomp_for_elvira = Comp::undefined();
      }
    }

    /* getter for normal vector method */
    inline NormMethod get_normal_vector_method_advance() {
      return norm_method_advance;
    }
    inline NormMethod get_normal_vector_method_curvature() {
      return norm_method_curvature;
    }

    /* setter for near-wall curvature method */
    void set_wall_curv_method(const CurvMethod wcm,
                              const Sign sig = Sign::undefined(),
                              const real cangle = -1.) {
      wall_curv_method = wcm;
      boil::oout<<"Wall curvature method: "<<wcm<<boil::endl;
      if(wcm==CurvMethod::HFmixedXZ()) {
        if       (sig==Sign::pos()) {
          mult_wall =  1;
        } else if(sig==Sign::neg()) {
          mult_wall = -1;
        } else {
          boil::oout<<"Phase at the origin must be specified! "
                    <<"Sign-neg: phi(orig)<0.5 and vice versa. "
                    <<"Exiting."
                    <<boil::endl;
          exit(0);
        }

        if(cangle>0.) {
          detachment_model.set_detachment_params(cangle);
        }
      }
      return;
    }

    /* getter for near-wall curvature method */
    inline CurvMethod get_wall_curv_method() { return wall_curv_method;};

    Vector fs;
    Vector * bndclr;
    Topology topo;

    Scalar nalpha;
    Scalar nx,ny,nz;/* normal to interface */
    Scalar adens;
    Scalar mx,my,mz;/* normal to interface, in real space */
  protected:
    virtual void ancillary(Scalar & scp);
    virtual void advance_x(Scalar & sca);
    virtual void advance_y(Scalar & sca);
    virtual void advance_z(Scalar & sca);

    void insert_bc_curv_divnorm();
    void insert_bc_curv_HFmixed(const Scalar & scp,
                                const Comp ctangential, const Comp cnormal,
                                const Sign sig);
    virtual real wall_curv_HFmixed_kernel(const real hc, const real hp,
                                          const real dc, const real dp,
                                          const real mult,
                                          const real cang);
    virtual real wall_curv_HFmixed_kernel(const real hm, const real hc, const real hp,
                                          const real dm, const real dc, const real dp,
                                          const real mult);
    
    void flood(Scalar & scp,const real mult);

    void cal_fs3(const Scalar & scp);
    void cal_fs_interp(const Scalar & scp);
    void curv_smooth();
    void extract_alpha(const Scalar & scp);
    void extract_alpha_near_bnd(const Scalar & scp);
    void fs_bnd(const Scalar & scp);
    void fs_bnd_nosubgrid(const Scalar & scp);
    void standardized_norm_vect(const Scalar & mx,
                                const Scalar & my,
                                const Scalar & mz,
                                Scalar & nx, Scalar & ny, Scalar & nz);
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void ib_norm(const Scalar & g);
    void ib_norm_cal(const int cc, const int i, const int j, const int k);
    void insert_bc(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_gradphi(const Scalar & g);
    void insert_bc_norm();
    void nib(const real & n1, const real & n2, const real & n3
           , const real & n4, const real & n5, const real & n6
           , real r[]);
    
    virtual void curv_HF();
    void curv_HF_kernel(arr3D & stencil, const arr3D & gridstencil,
                        const int imin, const int imax,
                        const real d1m, const real d1c, const real d1p,
                        const real d2m, const real d2c, const real d2p,
                        const real max_n,
                        const bool truedir1, const bool truedir2,
                        real & kap, int & flag);
                        //const int i, const int j, const int k);

    virtual void norm(const Scalar & color, const NormMethod & nm,
                      const bool extalp);

    void norm_cc(const Scalar & g);
    void norm_cc_kernel(real & nx_val, real & ny_val, real & nz_val,
                        real & nalpha_val,
                        Comp & mcomp,
                        const int i, const int j, const int k, 
                        const Scalar & sca);
    void norm_cc_kernel(real & nx_val, real & ny_val, real & nz_val, 
                        real & nalpha_val,
                        const int i, const int j, const int k, 
                        const Scalar & sca);
    void norm_cc_near_bnd(const Scalar & g);

    void norm_young(const Scalar & g);
    void norm_young_kernel(real & nx_val, real & ny_val, real & nz_val, 
                           real & nalpha_val,
                           const int i, const int j, const int k,
                           const Scalar & sca);

    void norm_mixed(const Scalar & g);
    void norm_mixed_kernel(real & nx_val, real & ny_val, real & nz_val, 
                           real & nalpha_val,
                           const int i, const int j, const int k,
                           const Scalar & sca);
  
    void bdnorm(Scalar & scp);
    void normal_vector_near_bnd(const Scalar & g, const NormMethod & nm);

    void normalize(real & r1, real & r2, real & r3);
    void normalize_l1(real & nx_l1, real & ny_l1, real & nz_l1,
                      const real nx, const real ny, const real nz);

    void nwall(const Scalar & g
             , const real & r1, const real & r2, const real & r3
             , const int & i1, const int & i2, const int &i3
             , real r[] );
    void smooth(const Scalar & sca, Scalar & scb, const int itnum);
    void true_norm_vect(const Scalar & nx,
                        const Scalar & ny,
                        const Scalar & nz,
                        Scalar & mx, Scalar & my, Scalar & mz);

    void update_at_walls(Scalar & scp);
    void wall_norm(const Scalar & sca);
    real kappa_ave(const real r1, const real r2);
    real kappa_ave(const real r1, const real r2, const int i1, const int i2);

    /* at the moment, these are NOT overwritten in the derived class */
    virtual real calc_v(const real r1, const real r2, const real r3, const real r4);
    virtual real calc_alpha(const real r1, const real r2, const real r3, const real r4);
    virtual real calc_flux(const real g, real c,
                           const real nx, const real ny, const real nz,
                           const real nalpha);

    void select_norm_cc(real & nx_val, real & ny_val, real & nz_val,
                        real & NxX, real & NyX, real & NzX,
                        real & NxY, real & NyY, real & NzY,
                        real & NxZ, real & NyZ, real & NzZ,
                        Comp * mcomp);
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
                       const Scalar & scp);

    void norm_cc_imin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_imax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmax(const Scalar &g, const int i,const int j, const int k);
    void vf_limiter();   

    /* elvira functions */
    void norm_elvira(const Scalar & sca, const bool extalp);
    void norm_elvira_kernel(real & nx_val, real & ny_val, real & nz_val,
                            real & nalpha_val,
                            const int i, const int j, const int k,
                            const Scalar & sca);
    void norm_elvira_kernel_full(real & nx_val, real & ny_val, real & nz_val,
                                 real & nalpha_val,
                                 const Comp & m, const int i, const int j, const int k,
                                 real valcc,real valmc,real valpc,real valcm,
                                 real valcp,real valmm,real valpm,real valmp,real valpp);
    void normalize_elvira(const real m, const real sig, 
                          real & nn1, real & nn2, real & nn3);
    real elvira_l2(const real alp,const real nn1,const real nn2,const real nn3,
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
    Scalar stmp;
    ScalarInt iflag,jflag;

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
    bool limit_color, use_subgrid, use_interp;
    real minclr, maxclr;

    Heaviside heavi;
    NormMethod norm_method_advance, norm_method_curvature;
    Comp mcomp_for_elvira;
    DetachmentModel detachment_model;
    CurvMethod bulk_curv_method, wall_curv_method;

    int nlayer;
    real cangle;
    real mult_wall;

    inline real signum(const real a, const real b) { return a*((b>0.)-(b<0.)); }
};	
#endif

