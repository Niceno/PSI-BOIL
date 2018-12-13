#include "phasechangevof.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::initialize() {

  phi=0.0;
  
  setflag();
 
  /* calculate normal vector */
  cal_norm_vect();

}

