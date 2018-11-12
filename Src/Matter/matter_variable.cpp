#include "matter.h"

/*============================================================================*/
void Matter::variable(const Set & s) {
/*--------------------------------------+
|  set certain property to be variable  |
+--------------------------------------*/

  if(s == Set::rho   ()) dens->varies(*dom);
  if(s == Set::mu    ()) visc->varies(*dom);
  if(s == Set::cp    ()) capa->varies(*dom);
  if(s == Set::lambda()) cond->varies(*dom);
  if(s == Set::gamma ()) diff->varies(*dom);
  if(s == Set::beta  ()) texp->varies(*dom);
  if(s == Set::sigma ()) 
    if(tens == NULL) {
      boil::oout << "# Fatal: using surface tension makes ";
      boil::oout << "sense only for mixtures. Exiting!"; 
      boil::oout << boil::endl;    
    } else
      tens->varies(*dom);
}

/*============================================================================*/
void Matter::variable() {
/*----------------------------------+
|  set all property to be variable  |
+----------------------------------*/

  dens->varies(*dom);
  visc->varies(*dom);
  capa->varies(*dom);
  cond->varies(*dom);
  diff->varies(*dom);
  texp->varies(*dom);
  if( tens!= NULL )
    tens->varies(*dom);
}
