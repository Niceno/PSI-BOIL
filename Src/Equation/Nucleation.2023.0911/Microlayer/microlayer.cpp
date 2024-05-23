#include "microlayer.h"
#include "../header.h"

/***************************************************************************//**
*  Constructor
*******************************************************************************/
Microlayer::Microlayer( Scalar & DM,
                        Scalar * MDOT,
                        Scalar * TPRS,
                        CommonHeatTransfer * CHT,
                        Heaviside * heavi,
                        const Times * t,
                        const real rs,
                        const real dmin, const real dmax,
                        Scalar * qsrc,
                        const Sign sig ) :
  Nucleation(CHT,heavi,t,rs,qsrc,sig),
  dmicro(&DM)
#ifndef USE_VOF_NUCL
  ,dSprev(*DM.domain())
#endif
{

  mdot = MDOT;
  tprs = TPRS;

  dmicro_min = dmin;
  dmicro_max = dmax;
  slope = 4.46e-3;  // Utaka's coefficient for water
  exp_slope = 1.0;
  rmax = boil::unreal;
  hresis = 0.0;//TIF::calculate_heat_transfer_resistance(cht->tifmodel.tref(), 
               //                                    rhov,mmass,latent,1.);
  boil::oout<<"microlayer:dmicro_min,dmicro_max= "<<dmin<<" "<<dmax<<"\n";

#ifndef USE_VOF_NUCL
  dmicro = boil::unreal;
  dSprev = dmicro.shape();
  str_dSprev = false;
#endif

  alloc1d ( & hflux_sum, 7);
  alloc1d ( & hflux_l, 7);
  alloc1d ( & hflux_v, 7);
  alloc1d ( & hflux_micro, 7);
  alloc1d ( & area_sum, 7);
  alloc1d ( & area_l, 7);
  alloc1d ( & area_v, 7);
  alloc1d ( & area_micro, 7);
}
