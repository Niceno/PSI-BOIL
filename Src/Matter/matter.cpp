#include "matter.h"

#define USE_PROSPERETTI

/*============================================================================*/
Matter::Matter(const Domain & d, const char * nm) {

  dom = & d;
  mixt = false;

  std::string sumnam;
  if(nm == NULL) {
    nam = "";
    sumnam = nam;
  } else {
    nam = nm;
    assert(nam.length() > 0);
    sumnam = nam + ".";
  }

  std::string prname; 
  prname = sumnam + "density";      dens = new Property(prname.c_str());
  prname = sumnam + "viscosity";    visc = new Property(prname.c_str());
  prname = sumnam + "capacity";     capa = new Property(prname.c_str());
  prname = sumnam + "conductivity"; cond = new Property(prname.c_str());
  prname = sumnam + "diffusivity";  diff = new Property(prname.c_str());
  prname = sumnam + "t_expansion";  texp = new Property(prname.c_str());
  prname = sumnam + "molar_mass";   molm = new Property(prname.c_str());
  tens = NULL;
  heat = NULL;

#ifdef USE_PROSPERETTI
  /* density over viscosity */
  dens_o_visc = new PropertyDiv(dens,visc);
  /* one over viscosity */
  //one_o_visc = new PropertyInv(visc);
#endif
}

/*============================================================================*/
Matter::Matter(const Matter & a, 
               const Matter & b, 
               const Scalar * ca,  
               const Scalar * cda,  
               const Scalar * cdb) {

  dom = ca->domain();

  mixt = true;

  assert(a.dens != NULL);
  assert(a.visc != NULL);
  assert(a.capa != NULL);
  assert(a.cond != NULL);
  assert(a.diff != NULL);
  assert(a.texp != NULL);
  assert(a.molm != NULL);
  assert(a.tens == NULL);
  assert(a.heat == NULL);
  assert(b.dens != NULL);
  assert(b.visc != NULL);
  assert(b.capa != NULL);
  assert(b.cond != NULL);
  assert(b.diff != NULL);
  assert(b.texp != NULL);
  assert(b.molm != NULL);
  assert(b.tens == NULL);
  assert(b.heat == NULL);
  if( ca->bc().count() == 0 ) {
    boil::oout << "# Fatal: defining mixture using concentration ";
    boil::oout << "variable without boundary conditions. Exiting!"; 
    boil::oout << boil::endl;    
    exit(0);
  }
  dens = new PropertyMix(a.dens, b.dens, ca, cda, cdb);
  capa = new PropertyMix(a.capa, b.capa, ca, cda, cdb);
  cond = new PropertyMix(a.cond, b.cond, ca, cda, cdb);
  diff = new PropertyMix(a.diff, b.diff, ca, cda, cdb);
  texp = new PropertyMix(a.texp, b.texp, ca, cda, cdb);
  molm = new PropertyMix(a.molm, b.molm, ca, cda, cdb);
  tens = new Property("surface-tension");
  heat = new Property("latent-heat");

#ifndef USE_PROSPERETTI
  visc = new PropertyMix(a.visc, b.visc, ca, cda, cdb);
#else /* force balance according to Prosperetti, 2002 */
  #if 1
  assert(a.dens_o_visc != NULL);
  assert(b.dens_o_visc != NULL);
  dens_o_visc = new PropertyMix(a.dens_o_visc,b.dens_o_visc,ca,cda,cdb);
  visc = new PropertyDiv(dens,dens_o_visc);
  #else
  assert(a.one_o_visc != NULL);
  assert(b.one_o_visc != NULL);
  one_o_visc = new PropertyMix(a.one_o_visc,b.one_o_visc,ca,cda,cdb);
  visc = new PropertyInv(one_o_visc);
  #endif
#endif

  if( a.nam.length() > 0 && b.nam.length() > 0 )
    nam = a.nam + "-" + b.nam;
  else
    nam = "mixture";
}

/*============================================================================*/
Matter::Matter(const Matter & a, 
               const Matter & b, 
               const Scalar * ca,  
               const Vector * bdca,  
               const Scalar * cda,  
               const Scalar * cdb) {

  dom = ca->domain();

  mixt = true;

  assert(a.dens != NULL);
  assert(a.visc != NULL);
  assert(a.capa != NULL);
  assert(a.cond != NULL);
  assert(a.diff != NULL);
  assert(a.texp != NULL);
  assert(a.molm != NULL);
  assert(a.tens == NULL);
  assert(a.heat == NULL);
  assert(b.dens != NULL);
  assert(b.visc != NULL);
  assert(b.capa != NULL);
  assert(b.cond != NULL);
  assert(b.diff != NULL);
  assert(b.texp != NULL);
  assert(b.molm != NULL);
  assert(b.tens == NULL);
  assert(b.heat == NULL);
  if( ca->bc().count() == 0 ) {
    boil::oout << "# Fatal: defining mixture using concentration ";
    boil::oout << "variable without boundary conditions. Exiting!"; 
    boil::oout << boil::endl;    
    exit(0);
  }
  dens = new PropertyMix(a.dens, b.dens, ca, bdca, cda, cdb);
  capa = new PropertyMix(a.capa, b.capa, ca, bdca, cda, cdb);
  cond = new PropertyMix(a.cond, b.cond, ca, bdca, cda, cdb);
  diff = new PropertyMix(a.diff, b.diff, ca, bdca, cda, cdb);
  texp = new PropertyMix(a.texp, b.texp, ca, bdca, cda, cdb);
  molm = new PropertyMix(a.molm, b.molm, ca, bdca, cda, cdb);
  tens = new Property("surface-tension");
  heat = new Property("latent-heat");

#if 1
  visc = new PropertyMix(a.visc, b.visc, ca, cda, cdb);
#else /* force balance according to Prosperetti, 2002 */
  assert(a.dens_o_visc != NULL);
  assert(b.dens_o_visc != NULL);
  dens_o_visc = new PropertyMix(a.dens_o_visc,b.dens_o_visc,ca,bdca,cda,cdb);
  visc = new PropertyDiv(dens,dens_o_visc);
#endif

  if( a.nam.length() > 0 && b.nam.length() > 0 )
    nam = a.nam + "-" + b.nam;
  else
    nam = "mixture";
}
