#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../Heaviside/MarchingSquares/MSaxisym/ms_axisym.h"
#include "../../Heaviside/MarchingCubes/marching_cubes.h"
#include "../../Topology/topology.h"
#include "../../../Global/global_realistic.h"
#include "../../../Parallel/communicator.h"
#include "../../../Ravioli/resrat.h"

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
        Krylov * S,
        Vector * bndclr = NULL);
    ~VOF();

    void new_time_step(){};
    void forward(Scalar & scp);

    void advance(const bool anci = true);
    void advance(Scalar & sca, const bool anci = true);

    void advance_with_extrapolation(const bool anci, const ResRat & resrat,
                                    const Vector & umixed, const Scalar & fext,
                                    const Matter * fluid_1, Vector * uvw_1,
                                    const Matter * fluid_2 = NULL,
                                    Vector * uvw_2 = NULL);

    void advance_with_extrapolation(Scalar & sca,
                                    const bool anci, const ResRat & resrat,
                                    const Vector & umixed, const Scalar & fext,
                                    const Matter * fluid_1, Vector * uvw_1,
                                    const Matter * fluid_2 = NULL,
                                    Vector * uvw_2 = NULL);

    void advance_phase_change(Scalar & scp);
    void advance_geometric(Scalar & scp);
    

    void curvature();
    void ancillary(); /* calcs ancillary params such as adens w/o advance */
    virtual void reconstruct_geometry();
    virtual void reconstruct_geometry(Scalar & scp);

    /* mainly used in VOFaxisym */
    virtual void color_to_vf(Scalar & color, Scalar & vf,
                             const bool extalp=true,const bool bdn=true) {
      vf = color;
      return;
    }

    void tension(Vector * vec, const Matter matt);
    void tension(Vector * vec, const Matter matt, const Scalar & scp);
    void totalvol();
    void front_minmax();
    void front_minmax( Range<real> xr
                     , Range<real> yr
                     , Range<real> zr );

    void init(){ ancillary(); };

    void extrapolate_velocity(const Scalar & scp, const Scalar & fext,
                              const Matter * fluid, const Vector & umixed, 
                              Vector & unew, const ResRat & resrat,
                              const Sign & sig, const bool flagging);

    virtual Scalar & color() {return phi;}
    const Vector & flow() { return vflow; }
    Heaviside * heaviside() { return heavi; }

    void interfacial_flagging(const Scalar & scp);

#include "vof_inline.h"

    Vector * bndclr;
    Topology topo;

    Scalar nalpha;
    Scalar nx,ny,nz;/* normal to interface */
  protected:
    void ancillary(Scalar & scp);
    virtual void advance_x(Scalar & sca);
    virtual void advance_y(Scalar & sca);
    virtual void advance_z(Scalar & sca);

    void ev_discretize(const Matter * fluid, const ScalarInt & pflag,
                       Matrix & A);
    void ev_flagging(const Scalar & scp, const ScalarInt & iflag,
                     ScalarInt & otpflag, const Sign & sig);
    bool ev_solve(const ScalarInt & pflag, const Matrix & A,
                  const Scalar & b, Scalar & x, Scalar & xold,
                  const bool init_guess, const int niter,
                  const ResRat & resrat);
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
    void insert_bc_curv_HFnormal(const Scalar & scp,
                                 const Comp ctangential, const Comp cnormal,
                                 const Sign sig, 
                                 const Range<int> ridx = Range<int>(-1,-2)); 
    void insert_bc_curv_HFparallel(const Scalar & scp,
                                   const Comp ctangential, const Comp cnormal,
                                   const Sign sig,
                                   const Range<int> ridx = Range<int>(-1,-2)); 
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
    
    void flood(Scalar & scp,const real mult);

    void cal_fs3(const Scalar & scp);
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
    Scalar stmp, stmp2;
    Scalar pold_pos, pold_neg; /* extrapolation */
    ScalarInt iflag,tempflag,tempflag2;

    Vector fs, vflow;
    Scalar adens;
    Scalar mx,my,mz;/* normal to interface, in real space */

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
    bool limit_color, use_interp, store_pressure_extrap;
    real minclr, maxclr;

    Heaviside * heavi;
    NormMethod norm_method_advance, norm_method_curvature;
    Comp mcomp_for_elvira;
    DetachmentModel detachment_model;
    CurvMethod bulk_curv_method, wall_curv_method;
    SubgridMethod subgrid_method;
    TopoMethod topo_method;
    HFset hf_set;

    int niter_pressure_extrap;
    int Nfilm_crit;
    real cangle;
    real mult_wall;

    inline real signum(const real a, const real b) { return a*((b>0.)-(b<0.)); }
    inline int signum(const int a, const int b) { return a*((b>0)-(b<0)); }
};	
#endif

