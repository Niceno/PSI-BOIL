#include "enthalpyfd.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
EnthalpyFD::EnthalpyFD(const Scalar & PHI, 
                       const Scalar & F,
                       const Vector & U,
                       const Vector & Uliq,
                       const Vector & Ugas,
                       Times & T,
                       Linear * S,
                       Matter * f,
                       const CommonHeatTransfer & CHT,
                       Matter * s) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S ),
  ftif   (  *PHI.domain()),
  cht(CHT),
  iflag(CHT.topo->iflag),
  iflagold(&(CHT.topo->iflagold)),
  uliq(&Uliq),
  ugas(&Ugas),
  flux_liq( *U.domain() ),
  flux_gas( *U.domain() ),
  bflag_struct(PHI),

  c_fff({ ConnectType::fluid,     ConnectType::fluid, ConnectType::fluid     }),
  c_sss({ ConnectType::solid,     ConnectType::solid, ConnectType::solid     }),

  c_iff({ ConnectType::interface, ConnectType::fluid, ConnectType::fluid     }),
  c_ffi({ ConnectType::fluid,     ConnectType::fluid, ConnectType::interface }),
  c_ifi({ ConnectType::interface, ConnectType::fluid, ConnectType::interface }),

  c_sff({ ConnectType::solid,     ConnectType::fluid, ConnectType::fluid     }),
  c_ffs({ ConnectType::fluid,     ConnectType::fluid, ConnectType::solid     }),
  c_sfs({ ConnectType::solid,     ConnectType::fluid, ConnectType::solid     }),

  c_sfi({ ConnectType::solid,     ConnectType::fluid, ConnectType::interface }),
  c_ifs({ ConnectType::interface, ConnectType::fluid, ConnectType::solid     }),

  c_ssf({ ConnectType::solid,     ConnectType::solid, ConnectType::fluid     }),
  c_fss({ ConnectType::fluid,     ConnectType::solid, ConnectType::solid     }),
  c_fsf({ ConnectType::fluid,     ConnectType::solid, ConnectType::fluid     }),

  c_ssi({ ConnectType::solid,     ConnectType::solid, ConnectType::interface }),
  c_iss({ ConnectType::interface, ConnectType::solid, ConnectType::solid     }),
  c_isi({ ConnectType::interface, ConnectType::solid, ConnectType::interface }),

  c_fsi({ ConnectType::fluid,     ConnectType::solid, ConnectType::interface }),
  c_isf({ ConnectType::interface, ConnectType::solid, ConnectType::fluid     })

{
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());
  laminar=true;
  ao_conv = AccuracyOrder::First();

  for_m(m) {
    flux_liq(m) = (*uliq)(m).shape();
    flux_gas(m) = (*ugas)(m).shape();
  }

  /* see header for explanation */
  if(solid()) {
    safe_solid = solid();
    if(cht.solid()==NULL) {
      boil::oout<<"enthalpyFD:Error!!! Solid is defined for EnthalpyFD ";
      boil::oout<<"but it is not defined for CommonHeat Transfer.\n";
      exit(0);
    }
  } else {
    safe_solid = fluid();
  }

  ftif = phi.shape();
  phi.bnd_update();

  convection_set(TimeScheme::forward_euler());
  diffusion_set(TimeScheme::backward_euler());

  discretize();
}
