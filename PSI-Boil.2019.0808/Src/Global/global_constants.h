#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <limits.h>
#include "global_precision.h"

#define DIM 3
 
namespace boil {
  /* big */
  extern real yotta;
  extern real zetta;
  extern real exa;
  extern real peta;
  extern real tera; 
  extern real giga; 
  extern real mega;
  extern real kilo;
  /* small */
  extern real milli;
  extern real micro;
  extern real nano; 
  extern real pico; 
  extern real femto;
  extern real atto;
  extern real zepto;
  extern real yocto;

  /* special */
  extern real unreal;
  extern int unint;

  /* pi */
  extern real pi;

  /* gravitational constant  */
  extern real g; 

  /* universal gas constant */
  extern real R;
}

#endif
