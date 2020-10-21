#ifndef CUSTOM_H
#define CUSTOM_H

#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Equation/Centered/CIPCSL2/cipcsl2.h"
#include "../Equation/Centered/VOF/vof.h"
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
                                    const Matter & sol,
                                    std::ostream & otp, const int NX);

  /* scalar error */
  real l1_scalar_error(const Scalar & sca, const Scalar & scb);
  real l2_scalar_error(const Scalar & sca, const Scalar & scb);
  real li_scalar_error(const Scalar & sca, const Scalar & scb);
  real l1_scalar_error_vol(const Scalar & sca, const Scalar & scb);
  real l2_scalar_error_vol(const Scalar & sca, const Scalar & scb);

  /* setup circle */
  void setup_circle_xz(Scalar & c, const real radius, const real xcent, const real zcent);

  /* setup sphere */
  real setup_sphere(Scalar & c, const real radius,
                    const real xcent, const real ycent, const real zcent,
                    const real mm = 20);

  /* setup square */
  void setup_plane(VOF & conc, const real nnx,
                   const real nny, const real nnz,
                   const real nalp);

  void setup_square_xz(VOF & conc, Scalar & tmp,
                       const real x0, const real z0,
                       const real sa, const real sb);

  /* prolongate plic */
  void prolongate_color_XZ(const VOF & concc, VOF & concf);
  inline real translate_v_coarse_to_fine(const int i_f, const int j_f, const int k_f,
                                         const int i_c, const int j_c, const int k_c,
                                         const Scalar & coarse, const Scalar & fine,
                                         const VOF & concc, VOF & concf);

  /* material properties as functions of temperature */
  inline real rho(const real rho0, const real beta, const real deltat) {
    return rho0/(1.+beta*deltat);
  } 

  template <class T>
  inline void print_line(T a) {
    boil::oout<<a<<boil::endl;
  }
}

#endif
