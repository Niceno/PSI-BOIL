#ifndef PHASECHANGE4_H
#define PHASECHANGE4_H

#include <cmath>
#include <set>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../CommonHeatTransfer/commonheattransfer.h"

#define IB

class PhaseChange4 : public Centered {
  public:
    PhaseChange4(const Scalar & mdot, 
                 const Scalar & mflx,
                 const Scalar & tprs,
                 const Scalar & vf,
                 const Scalar & vfs,
                 const Scalar & vs,
                 const Vector & u, 
                 CommonHeatTransfer & cht,
                 Times & t,
                 Matter * flu,
                 Matter * sol = NULL,
                 Sign matter_sig = Sign::pos());
 
    ~PhaseChange4();
    void update(const Scalar * diff_eddy = NULL);

    void mass_flux(const Scalar * diff_eddy = NULL);
    void initialize();
    void finalize();
    inline void sources() { sources_vfs(); sources_fext(); sources_sum(); }

#include "phasechange4_inline.h"

    void request_flux(const int i, const int j, const int k,
                      std::vector<real> & tv,
                      std::vector<real> & tl) const {
      tv = {txv[i][j][k],tyv[i][j][k],tzv[i][j][k],tnv[i][j][k]};
      tl = {txl[i][j][k],tyl[i][j][k],tzl[i][j][k],tnl[i][j][k]};
      return;
    }


  private:
    void m();
    real mdot_cut(const real mdotval, const real vfval,
                  const real rhol, real & mcut);
    real vfs_cut(real vfsval, real vfval);
    void mdot();

    void heat_flux(const Scalar * diff_eddy = NULL);
    void cal_hf(const Scalar * diff_eddy = NULL);
    void insert_bc_hf(const Scalar * diff_eddy);
    void dirac_source_terms();

    void sources_vfs();
    void sources_fext();
    void sources_sum();

    Scalar txv, tyv, tzv;
    Scalar txl, tyl, tzl;
    Scalar tnv, tnl;

    Scalar tprs, vf, vfs;
    Scalar nx;
    Scalar ny;
    Scalar nz;
    Scalar M;

    CommonHeatTransfer & cht;
    const Sign matter_sig; /* pos: liquid is phi=1 and vice versa */
 
    real smdot_pos, smdot_neg;
    
    AccuracyOrder accuracy_order;
    bool use_unconditional_extrapolation, discard_points_near_interface;
};	

#endif

