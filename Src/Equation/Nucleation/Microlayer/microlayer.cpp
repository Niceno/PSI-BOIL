#include "microlayer.h"

/***************************************************************************//**
*  Constructor
*******************************************************************************/
Microlayer::Microlayer( Scalar & DM,
                        Scalar * MDOT,
                        Scalar * VFS,
                        const Scalar * tpr,
                        Topology * topo,
                        Heaviside * heavi,
                        const TIF & TIFMODEL,
                        const Times * t,
                        Matter * f, const real rs, const real dm,
                        Matter * s,
                        Scalar * qsrc,
                        const Sign sig ) :
  Nucleation(topo,heavi,tpr,t,f,rs,qsrc,sig),
  dmicro(&DM),
  sld(s),
  dSprev(*DM.domain())
{

  mdot = MDOT;
  vfs = VFS;
  tifmodel = &TIFMODEL;

  dmicro_min = dm;
  dom = dmicro.domain();

  dSprev = dmicro.shape();
  
  slope = 4.46e-3;  // Utaka's coefficient for water
  exp_slope = 1.0;
  rmax = boil::unreal;

  hresis = Schrage::calculate_heat_transfer_resistance(tifmodel->tref(), 
                                                 rhov,mmass,latent);

  str_dSprev = false;

  alloc1d ( & hflux_micro, 7);
  alloc1d ( & area_sum, 7);
  alloc1d ( & area_l, 7);
  alloc1d ( & area_v, 7);

}
