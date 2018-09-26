#ifndef PHASECHANGE_H
#define PHASECHANGE_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Solver/Gauss/gauss.h"
#include "../../../Timer/timer.h"
#include "nucleation.h"

#define IB

class PhaseChange : public Centered {
  public:
    PhaseChange(const Scalar & s1, 
                const Scalar & s2,
                const Scalar & s3,
                const Scalar & s4,
                const Scalar & s5,
                const Scalar & s6,
                const Scalar & s7,
                const Vector & u, 
                Times & t,
                Matter * flu,
                real r1, real r2,
                Matter * sol = NULL,
                Nucleation * nucl = NULL);

 
    ~PhaseChange();
    void update(const Scalar * diff_eddy = NULL);
    void micro(Vector * vec, const Scalar * diff_eddy = NULL);
    void modify_vel(Vector & vec, const Scalar & sa,const Scalar & sb);

    // setter and getter for Mmicro and fmicro
    void set_Mmicro(real r){Mmicro=r;}
    void set_Fmicro(real r){Fmicro=r;}
    real get_Mmicro(){return Mmicro;}
    real get_Fmicro(){return Fmicro;}

    // setter and getter for interfacial thermal resistance
    void use_interfacial_resistance(bool b) {
      use_int_res = b;
      if (b) { boil::oout<<"Use interfacial thermal resistance"; 
               boil::oout<<"resint = "<<resint<<"\n";}
      else { boil::oout<<"Non-use interfacial thermal resistance"; }
    }
    void set_interfacial_resistance(real Tsat, real Rg, real alpha) {
      resint = Tsat * sqrt(2.0*pi*Rg*Tsat) / (latent*latent*rhov*alpha);
      boil::oout<<"PhaseChange:resint = "<<resint<<"\n";
      boil::oout<<"based on Tsat = "<<Tsat<<" (K), Rg = "<<Rg<<"\n";
      boil::oout<<"accomodation coefficient = "<<alpha<<"\n";
    }
    real get_interfacial_resistance() {
      return resint;
    }

    // heat flux
    real get_hflux(const Dir d){return hflux_total[int(d)];}
    real get_hflux_micro(const Dir d){return hflux_micro[int(d)];}
    real get_hflux_vapor(const Dir d){return hflux_vapor[int(d)];}
    real get_hflux_area(const Dir d){return area_sum[int(d)];}
    real get_hflux_area_l(const Dir d){return area_l[int(d)];}
    real get_hflux_area_v(const Dir d){return area_v[int(d)];}
    real get_smdot_pos(){return smdot_pos;}
    real get_smdot_neg(){return smdot_neg;}

    real get_turbP(){return turbP;}
    void set_turbP(real a){
      turbP=a;
      boil::oout<<"EnthalpyFD:turbP= "<<turbP<<"\n";
    }

  private:
    void cal_gradt(const Scalar * diff_eddy = NULL);
    void distfunc(const Scalar & sca, const int itnum);
    void dmicro_intermediate();
    bool incl_no_interface(const int i, const int j, const int k);
    //void init_dm0();
    void insert_bc_dist(Scalar & sca);
    void insert_bc_flag(ScalarInt & sca);
    void insert_bc_gradphic(const Scalar & sca);
    void insert_bc_norm();
    void gradphic(const Scalar & sca);
    void gradt (const Scalar * diff_eddy = NULL);
    void gradtx5( const int i, const int j, const int k,
              real * txm, real * txp, const int m) const;
    void gradty5( const int i, const int j, const int k,
              real * tym, real * typ, const int m) const;
    void gradtz5( const int i, const int j, const int k,
              real * tzm, real * tzp, const int m) const;

    void ext_gradt(Scalar & sca, const int iext);
    void m(const Scalar * diff_eddy = NULL);
    real marching_cube(const int i, const int j, const int k);
    real iso_length(const int i, const int j, const int k, const Dir d);
    void mdot_cut();
    void micro_shift();
    real mdot_cut(real m, real c);
    void mdot();
    void normalize(real &nx, real &ny, real &nz);
    real pcrate (real tpr);
    real pcfrc (real tpr);
    void qflux (const Scalar & sca);
    void setflag();
    void sources_clrs();
    void sources_fext();
    void sources_sum();
    void sources_tprs();
    void str_dSprev();
    void ib_set_dflag();
    void ib_ext_scalar(const Scalar & sca);
    void insert_bc(const Scalar & val);

    Nucleation * nucl;
    int nlayer,imodcal;
    Scalar tpr,tprs,clr,clrs,step;
    Scalar dist,stmp,stmp2,delta,M;
    ScalarInt dflag, iflag;
    Scalar nx,ny,nz;
    Scalar txv,tyv,tzv;
    Scalar txl,tyl,tzl;
    real dxmin,phisurf,pi,phimin,phimax;
    real latent, tsat, rhol, rhov, lambdal, lambdav, cpl, cpv;
    real rhoave, rhodlt;
    real epsl,epsnorm;
    real Mmicro,Fmicro;
    real smdot_pos, smdot_neg, smdot_pos_macro, smdot_neg_macro;
    real turbP;
    real * hflux_total, * hflux_micro, * hflux_vapor;
    real * area_sum, * area_l, * area_v;
    bool use_int_res;
    real resint;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: phasechange.h,v 1.15 2018/04/19 13:37:19 sato Exp $'/
+-----------------------------------------------------------------------------*/
