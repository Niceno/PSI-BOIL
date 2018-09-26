#include "levelset.h"
#include <cmath>

/******************************************************************************/
void LevelSet::advance() {
/***************************************************************************//**
*  \brief Compute convective equation by LevelSet method.
*******************************************************************************/
#ifdef DEBUG
  boil::oout<<"LevelSet_advance:start\n";
#endif

  convection();
  if(nredist>=1){
    if(time->current_step()%nredist==0) redist(16);
  }


#ifdef DEBUG
  boil::oout<<"LevelSet_advance:end\n";
#endif

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: levelset_advance.cpp,v 1.2 2012/09/13 08:42:26 niceno Exp $'/
+-----------------------------------------------------------------------------*/
