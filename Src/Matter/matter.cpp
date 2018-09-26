#include "matter.h"

/*============================================================================*/
Matter::Matter(const Domain & d) {

  dom = & d;

  nam = "";

  dens = new Property("density");
  visc = new Property("viscosity");
  capa = new Property("capacity");
  cond = new Property("conductivity");
  diff = new Property("diffusivity");
  texp = new Property("t_expansion");
  tens = NULL;
}

/*============================================================================*/
Matter::Matter(const Domain & d, const char * nm) {

  dom = & d;

  nam = nm;

  assert(nam.length() > 0);

  std::string prname; 
  prname = nam + "." + "density";      dens = new Property(prname.c_str());
  prname = nam + "." + "viscosity";    visc = new Property(prname.c_str());
  prname = nam + "." + "capacity";     capa = new Property(prname.c_str());
  prname = nam + "." + "conductivity"; cond = new Property(prname.c_str());
  prname = nam + "." + "diffusivity";  diff = new Property(prname.c_str());
  prname = nam + "." + "t_expansion";  texp = new Property(prname.c_str());
                                       tens = NULL;
}

/*============================================================================*/
Matter::Matter(const Matter & a, 
               const Matter & b, 
               const Scalar * ca,  
               const Scalar * cda,  
               const Scalar * cdb) {

  dom = ca->domain();

  assert(a.dens != NULL);
  assert(a.visc != NULL);
  assert(a.capa != NULL);
  assert(a.cond != NULL);
  assert(a.diff != NULL);
  assert(a.texp != NULL);
  assert(a.tens == NULL);
  assert(b.dens != NULL);
  assert(b.visc != NULL);
  assert(b.capa != NULL);
  assert(b.cond != NULL);
  assert(b.diff != NULL);
  assert(b.texp != NULL);
  assert(b.tens == NULL);
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
  tens = new Property("surface-tension");

  if( a.nam.length() > 0 && b.nam.length() > 0 )
    nam = a.nam + "-" + b.nam;
  else
    nam = "mixture";
}

/*-----------------------------------------------------------------------------+
 '$Id: matter.cpp,v 1.5 2015/08/19 12:01:43 badreddine Exp $'/
+-----------------------------------------------------------------------------*/
