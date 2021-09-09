#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "../../Parallel/mpi_macros.h"
#include "../../Global/global_realistic.h"
#include "../../Matter/matter.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Field/ScalarInt/scalarint.h"
#include "../../Field/Vector/vector.h"
#include "../../Domain/domain.h"
#include "../../Ravioli/accorder.h"
#include <set>

#include "topology_ravioli.h"

////////////////
//            //
//  Topology  //
//            //
////////////////
/* A class for safe argument passing from interface tracking to other classes */
class Topology {
  public:
    Topology(Scalar * VF, Scalar * CLR,
             Scalar * NX, Scalar * NY, Scalar * NZ, 
             Scalar * ADENS, Vector * FS, ScalarInt * IFLAG,
             const real clrsurf);
    
    ~Topology() {};

#include "topology_inline.h"

    /* store old variables */
    void new_time_step();

    void extrapolate(Scalar & sca, const Sign iext, const std::set<int> & testset);
    void extrapolate(Scalar & sca, const Sign iext, const std::set<int> & testset,
                     const ScalarInt & eflag);

    /* fs position */
     void cal_fs_interp(const Scalar & scp,Vector & fs,
                       const real tol_wall, const bool use_subgrid);

    void fs_bnd_nosubgrid(const Scalar & scp, Vector & fs,
                          const real & tol_wall);
    void fs_bnd_geometric(const Scalar & scp, Vector & fs,
                          const real & tol_wall);
    void fs_bnd_1D(const Scalar & scp, Vector & fs,
                   const real & tol_wall, const Sign & sig);
    
    /* capillary time step, coef is approx 1/sqrt(2pi) by default */
    static real capillary_ts(const Matter & mixed, const real dx,
                             const real coef = 0.3989);
    real wave_vel(const Matter & mixed,
                  const real dx, const real coef = 0.3989) const;
    real capillary_ts(const Matter & mixed, const Vector & vel,
                      const real coef = 0.3989) const;
    real capillary_ts(const Matter & mixed, 
                      const Vector & vel1, const Vector & vel2,
                      const real coef = 0.3989) const;

    /* interface boolean */
    bool interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k) const;
    bool interface_old(const Sign dir, const Comp m,
                       const int i, const int j, const int k) const;
    bool interface(const int i, const int j, const int k) const;

    /* distance to interface */
    real distance_int(const Sign dir, const Comp & m,
                      const int i, const int j, const int k,
                      Sign & cell_marker) const;

    real distance_int_x(const Sign dir,
                        const int i, const int j, const int k,
                        Sign & cell_marker) const;
    real distance_int_y(const Sign dir,
                        const int i, const int j, const int k,
                        Sign & cell_marker) const;
    real distance_int_z(const Sign dir,
                        const int i, const int j, const int k,
                        Sign & cell_marker) const;

    Sign distance1D_int_x(const int i, const int j, const int k,
                          const Sign dir, real & dist) const;
    Sign distance1D_int_y(const int i, const int j, const int k,
                          const Sign dir, real & dist) const;
    Sign distance1D_int_z(const int i, const int j, const int k,
                          const Sign dir, real & dist) const;

    /* old distance to interface */
    real distance_int_old(const Sign dir, const Comp & m,
                          const int i, const int j, const int k,
                          Sign & cell_marker) const;

    real distance_int_x_old(const Sign dir,
                            const int i, const int j, const int k,
                            Sign & cell_marker) const;
    real distance_int_y_old(const Sign dir,
                            const int i, const int j, const int k,
                            Sign & cell_marker) const;
    real distance_int_z_old(const Sign dir,
                            const int i, const int j, const int k,
                            Sign & cell_marker) const;

    Sign distance1D_int_x_old(const int i, const int j, const int k,
                              const Sign dir, real & dist) const;
    Sign distance1D_int_y_old(const int i, const int j, const int k,
                              const Sign dir, real & dist) const;
    Sign distance1D_int_z_old(const int i, const int j, const int k,
                              const Sign dir, real & dist) const;

    real zeroth_order_zeroth(const std::vector<StencilPoint> & stencil) const;
    real first_order_zeroth (const std::vector<StencilPoint> & stencil) const;
    real second_order_zeroth(const std::vector<StencilPoint> & stencil) const;
    real third_order_zeroth (const std::vector<StencilPoint> & stencil) const;
    real fourth_order_zeroth(const std::vector<StencilPoint> & stencil) const;
    real nth_order_zeroth(const std::vector<StencilPoint> & stencil,
                          const AccuracyOrder & order) const;

    real zeroth_order_first(const std::vector<StencilPoint> & stencil) const;
    real first_order_first (const std::vector<StencilPoint> & stencil) const;
    real second_order_first(const std::vector<StencilPoint> & stencil) const;
    real third_order_first (const std::vector<StencilPoint> & stencil) const;
    real fourth_order_first(const std::vector<StencilPoint> & stencil) const;
    real nth_order_first(const std::vector<StencilPoint> & stencil,
                         const AccuracyOrder & order) const;

    real zeroth_order_second(const std::vector<StencilPoint> & stencil) const;
    real first_order_second (const std::vector<StencilPoint> & stencil) const;
    real second_order_second(const std::vector<StencilPoint> & stencil) const;
    real third_order_second (const std::vector<StencilPoint> & stencil) const;
    real fourth_order_second(const std::vector<StencilPoint> & stencil) const;
    real nth_order_second(const std::vector<StencilPoint> & stencil,
                          const AccuracyOrder & order) const;

    real zeroth_order_third(const std::vector<StencilPoint> & stencil) const;
    real first_order_third (const std::vector<StencilPoint> & stencil) const;
    real second_order_third(const std::vector<StencilPoint> & stencil) const;
    real third_order_third (const std::vector<StencilPoint> & stencil) const;
    real fourth_order_third(const std::vector<StencilPoint> & stencil) const;
    real nth_order_third(const std::vector<StencilPoint> & stencil,
                         const AccuracyOrder & order) const;

    void zeroth_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil) const;
    void first_order_first_coefs (std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil) const;
    void second_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil) const;
    void third_order_first_coefs (std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil) const;
    void fourth_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil) const;
    void nth_order_first_coefs(std::vector<real> & coefs,
                               const std::vector<StencilPoint> & stencil,
                               const AccuracyOrder & order) const;

    /* testing */
    bool test_differences_first(const int count);
    bool test_differences_first(const std::vector<real> & stenpos,
                                const std::vector<real> & coefficients);
    real evaluate_polynomial(const int order,
                             const std::vector<real> & coefficients,
                             const real x);
    real evaluate_polynomial_derivative(const int order,
                             const std::vector<real> & coefficients,
                             const real x);

    /* front functions */
    void front_minmax(real * store_arr = NULL);
    void front_minmax(Range<real> xr,
                      Range<real> yr,
                      Range<real> zr,
                      real * store_arr = NULL);

    /* current variables */
    Scalar * vf;
    Scalar * clr;
    Scalar * nx;
    Scalar * ny;
    Scalar * nz;
    Scalar * adens;
    Vector * fs;
    ScalarInt * iflag;

    /* old variables */
    ScalarInt iflagold;
    Scalar clrold, vfold;
    Vector fsold;

  private:
    real frontPosition(const int i, const int j, const int k, const Comp m);

    int mmax_ext;
    real tol_ext, close_to_cc;
    real clrsurf;
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    
    ScalarInt stmp;
    Scalar delta, stmp2;

    const BndFlag bflag_struct;

    real dxmin;
};
#endif
