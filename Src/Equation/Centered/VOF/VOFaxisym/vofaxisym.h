#ifndef VOFAXISYM_H
#define VOFAXISYM_H

#include "../vof.h"

#define IB

using arr2D = std::vector< std::vector<real> >;

///////////
//       //
//  VOF  //
//       //
///////////
class VOFaxisym : public VOF {
  public:
    VOFaxisym(const Scalar & phi,
              const Scalar & f,
              const Scalar & kappa,
              const Vector & u, 
              Times & t,
              Krylov * S) :
              //Vector * bndclr = NULL) : /* bndclr not implemented! */
    VOF(phi,f,kappa,u,t,S,NULL),
    clr( *phi.domain() ),
    axistmp ( *phi.domain() ),
    Ktmp( *phi.domain() )
    {
       if(phi.domain()->is_cartesian()) {
         boil::oout<<"Warning: Initializing axisymmetric VOF on a Cartesian "
                   <<"domain!"<<boil::endl;
       }
       clr = phi.shape();
       axistmp = phi.shape();
       Ktmp = phi.shape();

       set_normal_vector_method_all(NormMethod::ElviraXZ());
       set_wall_curv_method(CurvMethod::HFmixedXZ(),Sign::pos());
       set_topo_method(TopoMethod::Hybrid());

       reconstruction_tolerance = 1e-4;
       reconstruction_maxiter = 5;

    }

    ~VOFaxisym() {};

    virtual void advance(const bool anci = true);
    virtual void advance(Scalar & sca, const bool anci = true);

    virtual void reconstruct_geometry();
    void color_to_vf(Scalar & color, Scalar & vf,
                     const bool extalp=true,const bool bdn=true);
    void vf_to_color(const Scalar & vf, Scalar & color);

    //real calc_alpha_axisymmetric(const real nnx, const real v, const real eta0);
    //real calc_v_axisymmetric(real nnx, real alp, real eta0, real & Kp);

    virtual Scalar & color() {return clr;}

    /* setter for reconstruction */
    void set_reconstruction_parameters(const real rtol, const int riter) {
      reconstruction_tolerance = rtol;
      reconstruction_maxiter = riter;
      boil::oout<<"set_reconstruction_parameters: tolerance= "<<rtol<<" ; "
                <<"maxiter= "<<riter<<"\n";
      return;
    }
    /* getter for reconstruction */
    inline real get_reconstruction_tolerance() { return(reconstruction_tolerance);}
    inline int get_reconstruction_maxiter() { return(reconstruction_maxiter);}

    real test_reconstruction(const Scalar & color, const Scalar & vf);

    //void forward_cartesian(Scalar & scp);
  protected:
    void forward_axisymmetric(const Scalar & color, Scalar & axip, Scalar & Kp);
    void backward_axisymmetric(const Scalar & vf, Scalar & alp);
    virtual void norm(const Scalar & color, const NormMethod & nm,
                      const bool extalp = true);

    virtual void advance_x(Scalar & sca);
    virtual void advance_y(Scalar & sca) {};
    virtual void advance_z(Scalar & sca);

    virtual void reconstruct_geometry(Scalar & scp);
    real linf_scalar_error(const Scalar & sca, const Scalar & scb
                          // ,int & i, int & j, int & k, bool & ebool
                           );
    
    virtual real wall_curv_HFmixed_kernel(const real hc, const real hp,
                                          const real dc, const real dp,
                                          const real mult,
                                          const real cang);
    virtual real wall_curv_HFmixed_kernel(const real hm, const real hc, const real hp,
                                          const real dm, const real dc, const real dp,
                                          const real mult);

    //virtual real calc_v(const real alpha, const real vma, const real vmb, const real vmc);
    real calc_v_axisymmetric(real nnx, real alp, real eta0, real & Kp);
    
    //virtual real calc_alpha(const real r1, const real r2,
    //                        const real r3, const real r4);
    real calc_alpha_axisymmetric(const real nnx, const real v, const real eta0);

    real calc_flux_axisymmetric(const real g,
                                const int i, const int j, const int k,
                                const Comp & mcomp);
    
    virtual void curv_HF();

#if 1
    void fill_stencil_x(arr2D & stencil, arr2D & gridstencil,
                        int & min, int & max,
                        const int i, const int j, const int k);
    void fill_stencil_z(arr2D & stencil, arr2D & gridstencil,
                        int & min, int & max,
                        const int i, const int j, const int k);

    void calculate_heights(arr2D & stencil, const arr2D & gridstencil,
                           const int imin, const int imax,
                           const real max_n, real & mult,
                           real & hm, real & hc, real & hp, real & nhc);

    void calculate_curvature_HF_axisymmetric(
                               const real hm, const real hc, const real hp,
                               const real dm, const real dc, const real dp,
                               const bool truedir, const real mult,
                               const Comp mcomp, const real xcent,
                               real & kap_cart, real & kap_cyl,
                               const int i, const int j, const int k);
#else
    void curv_HF_kernel_axisymmetric(
                               arr2D & stencil, const arr2D & gridstencil,
                               const int imin, const int imax,
                               const real dm, const real dc, const real dp,
                               const real max_n, const bool truedir,
                               const Comp mcomp, const real xcent,
                               real & kap, int & flag,
                               const int i, const int j, const int k
                                    );
#endif

    real xi_small_pos_triangle(real alpha, real mmx, real mmz);
    real xi_small_pos_trapezoid(real alpha, real mmx, real mmz);
    real xi_small_neg_triangle(real alpha, real mmx, real mmz);
    real xi_small_neg_trapezoid(real alpha, real mmx, real mmz);
    real xi_large_pos_triangle(real alpha, real mmx, real mmz);
    real xi_large_pos_trapezoid(real alpha, real mmx, real mmz);
    real xi_large_neg_triangle(real alpha, real mmx, real mmz);
    real xi_large_neg_trapezoid(real alpha, real mmx, real mmz);

    real phi_max_small_pos(const real mmx,const real mmz,const real eta0);
    real phi_max_small_neg(const real mmx,const real mmz,const real eta0);
    real phi_max_large_pos(const real mmx,const real mmz,const real eta0);
    real phi_max_large_neg(const real mmx,const real mmz,const real eta0);

    real phi_tr_small_pos(const real mmx,const real mmz,const real eta0);
    real phi_tr_small_neg(const real mmx,const real mmz,const real eta0);
    real phi_tr_large_pos(const real mmx,const real mmz,const real eta0);
    real phi_tr_large_neg(const real mmx,const real mmz,const real eta0);
    
    real phi_crit_pos(const real mmx,const real mmz,const real eta0);
    real phi_crit_neg(const real mmx,const real mmz,const real eta0);

    real alp_pos_triangle_subcrit(const real vf,const real mmx,const real mmz,
                                  const real eta0,const real phi_crit);
    real alp_pos_triangle_supercrit(const real vf,const real mmx,const real mmz,
                                    const real eta0,const real phi_crit);
    real alp_neg_triangle(const real vf,const real mmx,const real mmz,const real eta0);

    real alp_small_pos_trapezoid(const real vf,const real mmx,
                                 const real mmz,const real eta0);
    real alp_small_neg_trapezoid(const real vf,const real mmx,
                                 const real mmz,const real eta0);
    real alp_large_pos_trapezoid(const real vf,const real mmx,
                                 const real mmz,const real eta0);
    real alp_large_neg_trapezoid(const real vf,const real mmx,
                                 const real mmz,const real eta0);

    Scalar clr;
    Scalar axistmp, Ktmp;

    real reconstruction_tolerance;
    int reconstruction_maxiter;
};	
#endif

