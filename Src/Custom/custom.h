#ifndef CUSTOM_H
#define CUSTOM_H

#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Equation/Centered/CIPCSL2/cipcsl2.h"
#include "../Equation/Nucleation/nucleation.h"
#include "../Matter/matter.h"
#include <vector>
#include <iomanip>

////////////////////////
//                    //
//  Custom functions  //
//                    //
////////////////////////
namespace boil {

  /* cell-center velocities */
  void cell_center_velocities(const Vector & uvw,
                              Scalar & u, Scalar & v, Scalar & w);

  /* staggered velocities */
  void staggered_velocities(const Scalar & u, const Scalar & v, const Scalar & w,
                            Vector & uvw);
  /* backup io */
  std::string save_backup(const int ts, const bool irregular,
                   const Times & time,
                   const std::vector<Scalar*> & scalars,
                   const std::vector<std::string> & scalar_names,
                   const std::vector<Vector*> & vectors = {},
                   const std::vector<std::string> & vector_names = {},
                   const std::vector<Nucleation*> & nucls = {},
                   const std::vector<std::string> & nucl_names = {},
                   const std::vector<CIPCSL2*> & cipcsl2s = {},
                   const std::vector<std::string> & cipcsl2_names = {},
                   const std::vector<real*> & store_values = {});

  bool load_backup(const std::string & fname,
                   int & ts, Times & time,
                   const std::vector<Scalar*> & scalars,
                   const std::vector<std::string> & scalar_names,
                   const std::vector<Vector*> & vectors = {},
                   const std::vector<std::string> & vector_names = {},
                   const std::vector<Nucleation*> & nucls = {},
                   const std::vector<std::string> & nucl_names = {},
                   const std::vector<CIPCSL2*> & cipcsl2s = {},
                   const std::vector<std::string> & cipcsl2_names = {},
                   const std::vector<real*> & store_values = {});

  void rm_backup(const int ts,
                 const std::vector<Scalar*> & scalars,
                 const std::vector<std::string> & scalar_names,
                 const std::vector<Vector*> & vectors = {},
                 const std::vector<std::string> & vector_names = {},
                 const std::vector<Nucleation*> & nucls = {},
                 const std::vector<std::string> & nucl_names = {},
                 const std::vector<CIPCSL2*> & cipcsl2s = {},
                 const std::vector<std::string> & cipcsl2_names = {},
                 const std::vector<real*> & values = {});

  /* irun test and set */
  void test_irun(const std::string & testfile = "run.txt");
  void set_irun(const int val, const std::string & testfile = "run.txt");

  /* natural convection boundary layer thickness */
  real convective_boundary_layer_thickness(const Matter & mat,
                                           const real deltat);

  /* droplet parameters */
  void droplet_parameters_2D(const real cang, const real area,
                             real & radius, real & zcent, real & chord);
  void droplet_parameters_3D(const real cang, const real volume,
                             real & radius, real & zcent, real & chord);
  void droplet_parameters_3D(const real cang, real & volume,
                             real & radius, real & zcent, const real chord);
  void droplet_parameters_3D(const real cang, real & volume,
                             const real radius, real & zcent, real & chord);

  /* output scalar profile and wall heat transfer characteristics */
  void output_profile_xz(const Scalar & c, std::ostream & otp, 
                         const Range<int> RZ);
  void output_wall_heat_transfer_xz(const Scalar & tpr,
                                    const Vector & bndtpr,
                                    const real lambdas,
                                    std::ostream & otp, const int NX);

  /* scalar error */
  real l2_scalar_error(const Scalar & sca, const Scalar & scb);
  real li_scalar_error(const Scalar & sca, const Scalar & scb);

  /* setup circle */
  void setup_circle_xz(Scalar & c, const real radius, const real xcent, const real zcent);

  /* setup sphere */
  real setup_sphere(Scalar & c, const real radius,
                    const real xcent, const real ycent, const real zcent,
                    const real mm = 20);

  /* interpolate velocities */
  void interpolateXZ(const Vector & coarse, Vector & fine, 
                     const Scalar & cs, const Scalar & mdot,
                     const Matter & fluid);

  void interpolate_pressure_solve_2D(
                     real & x, real & y, real & z, real & q,
                     const real & F_x, const real & F_y,
                     const real & F_z, const real & F_q,
                     const real & c_a, const real & c_b,
                     const real & c_c, const real & c_d,
                     const real & Q_a, const real & Q_b,
                     const real & Q_c, const real & Q_d,
                     const real & Q_wb, const real & Q_wt,
                     const real & Q_eb, const real & Q_et,
                     const real & Q_bw, const real & Q_be,
                     const real & Q_tw, const real & Q_te);

  /* restrict scalar */
  void restrictXZ(const Scalar & fine, Scalar & coarse);
  void restrictXZ_area(const Scalar & fine, Scalar & coarse);

  /* restrict vector */
  void restrictXZ_vector_simple(const Vector & fine, Vector & coarse,
                                const Scalar & fs, const Scalar & cs);

  /* material properties as functions of temperature */
  inline real rho(const real rho0, const real beta, const real deltat) {
    return rho0/(1.+beta*deltat);
  } 
}

#endif
