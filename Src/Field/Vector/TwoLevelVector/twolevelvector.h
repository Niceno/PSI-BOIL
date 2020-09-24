#ifndef TWOLEVELVECTOR_H
#define TWOLEVELVECTOR_H

#include "../vector.h"
#include "../../../Matter/matter.h"

////////////////////////
//                    //
//  Two-level Vector  //
//                    //
////////////////////////
class TwoLevelVector {
  public:
    explicit TwoLevelVector(const TwoLevelDomain & d) :
      coarse(d.coarse()),
      fine(d.fine()),
      levels({&coarse,&fine})
    {
    }

    void restrict_area_XZ(const Scalar & cs);
    void divergence_free_interpolate_XZ(
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

    Vector coarse;
    Vector fine;

    Vector * operator [] (const int l) {return levels[l];}

    std::array<Vector*,2> levels;
};
#endif
