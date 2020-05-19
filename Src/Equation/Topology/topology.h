#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "../../Parallel/mpi_macros.h"
#include "../../Global/global_realistic.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Field/ScalarInt/scalarint.h"
#include "../../Field/Vector/vector.h"
#include "../../Domain/domain.h"
#include <set>

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

    /* store old variables */
    void new_time_step();

    void extrapolate(Scalar & sca, const Sign iext, const std::set<int> & testset);
    void extrapolate(Scalar & sca, const Sign iext, const std::set<int> & testset,
                     const ScalarInt & eflag);
    
    inline int get_extrapolation_iters() const { return mmax_ext; }
    inline real get_extrapolation_tol() const  { return tol_ext; }
    inline void set_extrapolation_params(const int mnew, const real tolnew) {
      mmax_ext = mnew;
      tol_ext = tolnew;
      boil::oout<<"Topology::extrapolationparams: "<<mnew<<" "<<tolnew<<"\n";
    }

    /* interface boolean */
    bool interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);
    bool interface_old(const Sign dir, const Comp m,
                       const int i, const int j, const int k);
    bool interface(const int i, const int j, const int k);

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

    inline bool above_interface(const int i, const int j, const int k)
      const { return (*clr)[i][j][k]>=clrsurf; }
    inline bool above_interface_old(const int i, const int j, const int k)
      const { return clrold[i][j][k]>=clrsurf; }
    inline bool under_interface(const int i, const int j, const int k)
      const { return (*clr)[i][j][k]<clrsurf; }
    inline bool under_interface_old(const int i, const int j, const int k)
      const { return clrold[i][j][k]<clrsurf; }

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
    int mmax_ext;
    real tol_ext, close_to_cc;
    real clrsurf;
    
    ScalarInt stmp;
    Scalar delta, stmp2;
};
#endif
