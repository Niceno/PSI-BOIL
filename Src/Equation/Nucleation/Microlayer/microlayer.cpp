#include "microlayer.h"

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
  dmicro(&DM)//,
  //dSprev(*DM.domain())
{

  mdot = MDOT;
  tprs = TPRS;

  dmicro_min = dmin;
  dmicro_max = dmax;

  //dSprev = dmicro.shape();
  
  slope = 4.46e-3;  // Utaka's coefficient for water
  exp_slope = 1.0;
  rmax = boil::unreal;

  hresis = TIF::calculate_heat_transfer_resistance(cht->tifmodel.tref(), 
                                                   rhov,mmass,latent,1.);

  //str_dSprev = false;

  alloc1d ( & hflux_micro, 7);
  alloc1d ( & area_sum, 7);
  alloc1d ( & area_l, 7);
  alloc1d ( & area_v, 7);

}
