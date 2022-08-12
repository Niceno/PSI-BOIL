#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../Heaviside/MarchingSquares/MSaxisym/ms_axisym.h"
#include "../../Heaviside/MarchingCubes/marching_cubes.h"
#include "../../Topology/topology.h"
#include "../../../Global/global_realistic.h"
#include "../../../Global/global_func.h"
#include "../../../Parallel/communicator.h"
#include "../../../Ravioli/restol.h"

#define IB

using arr2D = std::vector< std::vector<real> >;
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
        Linear * S,
        Vector * bndclr = NULL);
    ~VOF();

    void new_time_step(const Scalar * diff_eddy = NULL) {
      if(!is_initialized)
        init();
      topo->new_time_step(); 
    };
    void forward(Scalar & scp);

    void advance(const bool anci = true);
    void advance(const Scalar & sca, const bool anci = true);

    void advance_with_extrapolation(const bool anci, const ResTol & restol,
                                    const Vector & umixed, const Scalar & fext,
                                    const Matter * fluid_1, Vector * uvw_1,
                                    const Matter * fluid_2 = NULL,
                                    Vector * uvw_2 = NULL);

    void advance_with_extrapolation(const Scalar & sca,
                                    const bool anci, const ResTol & restol,
                                    const Vector & umixed, const Scalar & fext,
                                    const Matter * fluid_1, Vector * uvw_1,
                                    const Matter * fluid_2 = NULL,
                                    Vector * uvw_2 = NULL);

    void advance_with_extrapolation(const bool anci, const ResTol & restol,
                                    const Vector & umixed,
                                    const Matter * fluid_1, Vector * uvw_1,
                                    const Matter * fluid_2 = NULL,
                                    Vector * uvw_2 = NULL);

    void advance_with_extrapolation(const Scalar & sca,
                                    const bool anci, const ResTol & restol,
                                    const Vector & umixed,
                                    const Matter * fluid_1, Vector * uvw_1,
                                    const Matter * fluid_2 = NULL,
                                    Vector * uvw_2 = NULL);

    void advance_phase_change(Scalar & scp);
    void advance_geometric(Scalar & scp);
    

    void curvature();
    /* calcs ancillary params such as adens w/o advance */
    void ancillary(const bool reconstruct = true);
    virtual void reconstruct_geometry();
    virtual void reconstruct_geometry(Scalar & scp);

    /* mainly used in VOFaxisym */
    virtual void color_to_vf(Scalar & color, Scalar & vf,
                             const bool nvec=true,const bool extalp=true,
                             const bool bdn=true) {
      vf = color;
      return;
    }

    virtual void color_to_vf(const bool nvec=true,const bool extalp=true,
                             const bool bdn=true) {
      return;
    }

    void tension(Vector * vec, const Matter & matt);
    void tension(Vector * vec, const Matter & matt, const Scalar & scp);
    real totalvol(real * vaps = NULL);
    real totalvol(Range<real> xr
                , Range<real> yr
                , Range<real> zr
                , real * vaps = NULL);
    void front_minmax(const std::string & nm = "x-min-front=");
    void front_minmax( Range<real> xr
                     , Range<real> yr
                     , Range<real> zr
                     , const std::string & nm = "x-min-front=");

    void init(){ ancillary(); is_initialized = true; };

    void extrapolate_velocity(const Scalar & scp, const Scalar & fext,
                              const Matter * fluid, const Vector & umixed, 
                              Vector & unew, const ResTol & restol,
                              const Sign & sig, const bool flagging);
    void extrapolate_velocity(const Scalar & scp,
                              const Matter * fluid, const Vector & umixed,
                              Vector & unew, const ResTol & restol,
                              const Sign & sig, const bool flagging);

    virtual Scalar & color() {return phi;}
    virtual const Scalar & color() const {return phi;}
    const Vector & flow() { return vflow; }
    Heaviside * heaviside() { return heavi; }

    void interfacial_flagging(const Scalar & scp);

    real output_cangle_2d(const Comp ctangential, const Comp cnormal,
                          const Sign sig, const Sign mult_wall);
   
    real extract_cl_velocity_2d(const Comp ctangential, const Comp cnormal,
                                const Sign sig,
                                int * IG = NULL, int * PN = NULL,
                                const Range<int> ridx = Range<int>(-1,-2));

    real translate_v(const int i, const int j, const int k,
                     const real dx, const real dy, const real dz,
                     const real fx, const real fy, const real fz,
                     real & nnx, real & nny, real & nnz, real & naa,
                     const Scalar & scp) const;

    /* only for comparison purposes!!! */
    void cal_adens_geom(Scalar & adensgeom, const Scalar & sca,
                        const bool use_vicinity = true);
    void cal_adens_gradclr(Scalar & adensgeom, const Scalar & sca);
    void cal_adens_gradclr_2phi(Scalar & adensgeom, const Scalar & sca);
    void cal_adens_gradclr_6phi(Scalar & adensgeom, const Scalar & sca);

    void color_minmax();
#include "vof_inline.h"

    Vector * bndclr;
    Topology * topo;
    DetachmentModel detachment_model;

    Scalar nalpha;
    Scalar nx,ny,nz;/* normal to interface */
  protected:
    void ancillary(Scalar & scp, const bool reconstruct = true);
    virtual void advance_x(const Scalar & sca, Scalar & cellvol);
    virtual void advance_y(const Scalar & sca, Scalar & cellvol);
    virtual void advance_z(const Scalar & sca, Scalar & cellvol);

    void update_phi(const Scalar & cellvol, Scalar & sca);
    void advect_naive(Scalar & scp);
    void advect_reconstructed(Scalar & scp);
    void advect_bounded(Scalar & scp);
    void divergence_skew(const Comp & m, const ScalarInt & marker,
                         Scalar & cellvol);

    void ev_discretize(const Matter * fluid, const ScalarInt & pflag,
                       Matrix & A);
    void ev_flagging(const Scalar & scp, const ScalarInt & iflag,
                     ScalarInt & otpflag, const Sign & sig);
    bool ev_solve(const ScalarInt & pflag, const Matrix & A,
                  const Scalar & b, Scalar & x, Scalar & xold,
                  const bool init_guess, const int niter,
                  const ResTol & restol);
    void ev_calculate_source(const ScalarInt & pflag, const Vector & u,
                             Scalar & fext);
    void ev_project(const ScalarInt & pflag, const Matter * fluid,
                    const Scalar & frc, Vector & u);
    void ev_complement(const ScalarInt & pflag, const Scalar & scp,
                       const Vector & umixed, const Vector & uliq, 
                       Vector & ugas, const Vector * bndclr);

    //void interfacial_flagging(const Scalar & scp);
    bool Interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);
    bool Interface(const int i, const int j, const int k);

    void insert_bc_curv_divnorm();
    void wall_norm(const Scalar & sca);
    void nwall(const Scalar & g,
               const real & r1, const real & r2, const real & r3,
               const int & i1, const int & i2, const int &i3,
               real r[] );

#if 0
    /* 2D-variants */
    void insert_bc_curv_HFnormal(const Scalar & scp,
                                 const Comp ctangential, const Comp cnormal,
                                 const Sign sig);
                                // const Range<int> ridx = Range<int>(-1,-2)); 
    void insert_bc_curv_HFparallel(const Scalar & scp,
                                   const Comp ctangential, const Comp cnormal,
                                   const Sign sig);
                                   //const Range<int> ridx = Range<int>(-1,-2)); 
    void insert_bc_curv_HFmixed(const Scalar & scp,
                                const Comp ctangential, const Comp cnormal,
                                const Sign sig);

    virtual real wall_curv_HFnormal_kernel(const real x0, const real hm,
                                           const real hc, const real hp,
                                           const real dm,
                                           const real dc, const real dp,
                                           const real mult, const real cang);
    virtual real wall_curv_HFparallel_kernel(const real hc, const real hp,
                                             const real dc, const real dp,
                                             const real mult, const real cang);
    virtual real wall_curv_HFparallel_kernel(const real hm,
                                             const real hc, const real hp,
                                             const real dm,
                                             const real dc, const real dp,
                                             const real mult);

    /* 3D-variant */
    void insert_bc_curv_HFnormal(const Scalar & scp,
                                 const Comp cnormal,
                                 const Sign sig);

    real wall_curv_HFnormal_kernel(arr2D & heights,
                                   const arr2D & distances,
                                   const std::vector<real> & dx,
                                   const std::vector<real> & dy,
                                   const bool truedir1, const bool truedir2,
                                   const real mult, const real max_n,
                                   const real cang);
#endif
    
    void flood(Scalar & scp,const real mult);
    void output_cangle_2d(const Scalar & scp,
                          const Comp ctangential, const Comp cnormal,
                          const Sign sig, const Sign mult_wall,
                          const Range<int> ridx,
                          real & h0, real & h1, real & h2,
                          real & dzzt0, real & dzzc0,
                          real & dzzt1, real & dzzc1);

    void cal_fs3(const Scalar & scp);
    void fs_bnd_symmetry(const Scalar & scp, Vector & fs,
                         const real & tol_wall);
    void extract_alpha(const Scalar & scp);
    void extract_alpha_near_bnd(const Scalar & scp);
    void standardized_norm_vect(const Scalar & mx,
                                const Scalar & my,
                                const Scalar & mz,
                                Scalar & nx, Scalar & ny, Scalar & nz);
    void gradphi(const Scalar & g);
    void ib_norm(const Scalar & g);
    void ib_norm_cal(const int cc, const int i, const int j, const int k);
    void insert_bc_gradphi(const Scalar & g);
    void nib(const real & n1, const real & n2, const real & n3
           , const real & n4, const real & n5, const real & n6
           , real r[]);
    
    void curv_HF();
    void curv_HF_kernel(arr3D & stencil, const arr3D & gridstencil,
                        const arr2D & wall_indicator, const Comp & mcomp,
                        const int imin, const int imax,
                        const real d1m, const real d1c, const real d1p,
                        const real d2m, const real d2c, const real d2p,
                        const real cang_m, const real cang_p,
                        const real max_n, const real nnj, const real nnk,
                        const bool truedir1, const bool truedir2,
                        const real xcent,
                        real & kap, int & flag
                        ,const int i, const int j, const int k
                       ) const;
    virtual real calculate_curvature_HF(const arr2D & heights,
                                const real d1m, const real d1c, const real d1p,
                                const real d2m, const real d2c, const real d2p,
                                const bool truedir1, const bool truedir2,
                                const real mult, const real max_n,
                                const real xcent) const;

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

    void smooth(Scalar & sca, const int itnum);
    void true_norm_vect(const Scalar & nx,
                        const Scalar & ny,
                        const Scalar & nz,
                        Scalar & mx, Scalar & my, Scalar & mz);

    void update_at_walls(Scalar & scp);
    void update_at_walls_custom(Scalar & scp);
    real kappa_ave(const real r1, const real r2);
    real kappa_ave(const real r1, const real r2, const int i1, const int i2);

    /* at the moment, these are NOT overwritten in the derived class */
    virtual real calc_v(const real r1, const real r2, const real r3, const real r4) const ;
    virtual real calc_alpha(const real r1, const real r2, const real r3, const real r4) const;
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

    Scalar kappa;        /* curvature */
    Scalar stmp, stmp2;
    Scalar pold_pos, pold_neg; /* extrapolation */
    ScalarInt iflag,tempflag,tempflag2;

    Vector fs, vflow;
    Scalar adens;
    Scalar mx,my,mz;/* normal to interface, in real space */

    Matter jelly;   /* virtual fluid needed by parent constructor */
    real theta;
    real epsnorm;
    real phisurf;
    real tol_wall, tol_flux, tol_ext, flux_cfl;
    real ww, dxmin;
    const BndFlag bflag_struct;
    bool limit_color, use_interp, store_pressure_extrap;
    bool is_initialized;
    real minclr, maxclr;

    Heaviside * heavi;
    NormMethod norm_method_advance, norm_method_curvature;
    Comp mcomp_for_elvira, wall_curv_dir;
    CurvMethod bulk_curv_method, wall_curv_method;
    SubgridMethod subgrid_method;
    AdvectionMethod advect_method;
    HFset hf_set;

    boil::func_ijk_real cangle_func;
    boil::func_scijk_real update_at_walls_func;

    /* labels for advection rotation */
    std::array<int,6> label_adv;

    int niter_pressure_extrap;
    real cangle0;

    bool cangle_variable, update_at_walls_variable, extrapolate_ib;
};	
#endif

