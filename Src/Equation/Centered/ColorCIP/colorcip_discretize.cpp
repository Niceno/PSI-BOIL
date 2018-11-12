#include "colorcip.h"

/******************************************************************************/
void ColorCIP::discretize() {

  /* correct on the boundaries */
  create_system_bnd();   

}
