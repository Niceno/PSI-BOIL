#include "microlayer.h"

/***************************************************************************//**
*  Constructor
*******************************************************************************/
Microlayer::Microlayer( Scalar & DM,
                        Scalar * MDOT,
                        Scalar * TPRS,
                        const Scalar * tpr,
                        Topology * topo,
                        Heaviside * heavi,
                        const TIF & TIFMODEL,
                        const Times * t,
                        Matter * f, const real rs,
                        const real dmin, const real dmax,
                        Matter * s,
                        Scalar * qsrc,
                        const Sign sig ) :
  Nucleation(topo,heavi,tpr,t,f,rs,qsrc,sig),
  dmicro(&DM),
  sld(s),
  dSprev(*DM.domain())
{

  mdot = MDOT;
  tprs = TPRS;
  tifmodel = &TIFMODEL;

  dmicro_min = dmin;
  dmicro_max = dmax;
  dom = dmicro.domain();

  dSprev = dmicro.shape();
  
  slope = 4.46e-3;  // Utaka's coefficient for water
  exp_slope = 1.0;
  rmax = boil::unreal;

  hresis = Schrage::calculate_heat_transfer_resistance(tifmodel->tref(), 
                                                 rhov,mmass,latent,1.);

  str_dSprev = false;

  alloc1d ( & hflux_micro, 7);
  alloc1d ( & area_sum, 7);
  alloc1d ( & area_l, 7);
  alloc1d ( & area_v, 7);

}
