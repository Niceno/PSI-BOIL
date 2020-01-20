#include "matter.h"

/*============================================================================*/
Matter::Matter(const Domain & d) {

  dom = & d;

  nam = "";

  mixt = false;

  dens = new Property("density");
  visc = new Property("viscosity");
  capa = new Property("capacity");
  cond = new Property("conductivity");
  diff = new Property("diffusivity");
  texp = new Property("t_expansion");
  molm = new Property("molar_mass");
  tens = NULL;
  heat = NULL;
}

/*============================================================================*/
Matter::Matter(const Domain & d, const char * nm) {

  dom = & d;

  nam = nm;

  mixt = false;

  assert(nam.length() > 0);

  std::string prname; 
  prname = nam + "." + "density";      dens = new Property(prname.c_str());
  prname = nam + "." + "viscosity";    visc = new Property(prname.c_str());
  prname = nam + "." + "capacity";     capa = new Property(prname.c_str());
  prname = nam + "." + "conductivity"; cond = new Property(prname.c_str());
  prname = nam + "." + "diffusivity";  diff = new Property(prname.c_str());
  prname = nam + "." + "t_expansion";  texp = new Property(prname.c_str());
  prname = nam + "." + "molar_mass";   molm = new Property(prname.c_str());
                                       tens = NULL;
                                       heat = NULL;
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
  visc = new PropertyMix(a.visc, b.visc, ca, cda, cdb);
  capa = new PropertyMix(a.capa, b.capa, ca, cda, cdb);
  cond = new PropertyMix(a.cond, b.cond, ca, cda, cdb);
  diff = new PropertyMix(a.diff, b.diff, ca, cda, cdb);
  texp = new PropertyMix(a.texp, b.texp, ca, cda, cdb);
  molm = new PropertyMix(a.molm, b.molm, ca, cda, cdb);
  tens = new Property("surface-tension");
  heat = new Property("latent-heat");

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
  visc = new PropertyMix(a.visc, b.visc, ca, bdca, cda, cdb);
  capa = new PropertyMix(a.capa, b.capa, ca, bdca, cda, cdb);
  cond = new PropertyMix(a.cond, b.cond, ca, bdca, cda, cdb);
  diff = new PropertyMix(a.diff, b.diff, ca, bdca, cda, cdb);
  texp = new PropertyMix(a.texp, b.texp, ca, bdca, cda, cdb);
  molm = new PropertyMix(a.molm, b.molm, ca, bdca, cda, cdb);
  tens = new Property("surface-tension");
  heat = new Property("latent-heat");

  if( a.nam.length() > 0 && b.nam.length() > 0 )
    nam = a.nam + "-" + b.nam;
  else
    nam = "mixture";
}
