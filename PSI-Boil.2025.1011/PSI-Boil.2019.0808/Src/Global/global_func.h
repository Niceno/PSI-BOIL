#ifndef FUNC_H
#define FUNC_H

#include <functional>
#include "global_precision.h"
#include "../Ravioli/comp.h"
#include "../Ravioli/sign.h"

namespace boil {
  typedef std::function<real(const int, const int, const int)> func_ijk_real;
  typedef std::function<real(const Sign, const Comp, const int, const int, const int)> func_scijk_real;
}

#endif
