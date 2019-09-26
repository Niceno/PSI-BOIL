#ifndef VOFAXISYM_H
#define VOFAXISYM_H

#include "../vof.h"

#define IB

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
              Krylov * S,
              Vector * bndclr = NULL) :
    VOF(phi,f,kappa,u,t,S,bndclr),
    clr( *phi.domain() ),
    axistmp ( *phi.domain() ),
    axistmp2( *phi.domain() ),
    Ktmp( *phi.domain() )
    {
       if(phi.domain()->is_cartesian()) {
         boil::oout<<"Warning: Initializing axisymmetric VOF on a Cartesian "
                   <<"domain!"<<boil::endl;
       }
       clr = phi.shape();
       axistmp = phi.shape();
       axistmp2= phi.shape();
       Ktmp = phi.shape();

       set_normal_vector_method_all(NormMethod::ElviraXZ());

       reconstruction_tolerance = 1e-4;
       reconstruction_maxiter = 20;
    }

    ~VOFaxisym() {};

    virtual void advance(const bool anci = true);
    virtual void advance(Scalar & sca, const bool anci = true);

    virtual void ancillary() {};
    void reconstruct_geometry();
    void color_to_vf(const Scalar & color, Scalar & vf);
    void vf_to_color(const Scalar & vf, Scalar & color);

    //real calc_alpha_axisymmetric(const real nnx, const real v, const real eta0);
    //real calc_v_axisymmetric(real nnx, real alp, real eta0, real & Kp);

    Scalar & color() {return clr;}

    /* setter for reconstruction */
    void set_reconstruction_parameters(const real rtol, const int riter) {
      reconstruction_tolerance = rtol;
      reconstruction_maxiter = riter;
      boil::oout<<"set_reconstruction_parameters: tolerance= "<<rtol<<" ; "
                <<"maxiter= "<<riter<<"\n";
      return;
    }
    /* getter for reconstruction */
    real get_reconstruction_tolerance() { return(reconstruction_tolerance);}
    int get_reconstruction_maxiter() { return(reconstruction_maxiter);}

    real test_reconstruction(const Scalar & color, const Scalar & vf);

  protected:
    void forward_axisymmetric(const Scalar & color, Scalar & axip, Scalar & Kp);
    void backward_axisymmetric(const Scalar & vf, Scalar & alp);
    void norm_axisymmetric(const Scalar & color);

    virtual void advance_x(Scalar & sca);
    virtual void advance_y(Scalar & sca) {};
    virtual void advance_z(Scalar & sca);

    virtual void ancillary(Scalar & scp) {};
    void reconstruct_geometry(const Scalar & scp);
    real linf_scalar_error(const Scalar & sca, const Scalar & scb);

    //virtual real calc_v(real alpha, real vma, real vmb, real vmc);
    real calc_v_axisymmetric(real nnx, real alp, real eta0, real & Kp);
    
    //virtual real calc_alpha(const real r1, const real r2,
    //                        const real r3, const real r4);
    real calc_alpha_axisymmetric(const real nnx, const real v, const real eta0);

    real calc_flux_axisymmetric(const real g,
                                const int i, const int j, const int k,
                                const Comp & mcomp);

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
    Scalar axistmp, axistmp2, Ktmp;

    real reconstruction_tolerance;
    int reconstruction_maxiter;
};	
#endif

