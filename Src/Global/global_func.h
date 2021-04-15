#ifndef FUNC_H
#define FUNC_H

#include <functional>
#include "global_precision.h"

namespace boil {
  typedef std::function<real(const int, const int, const int)> func_ijk_real;
}

#endif
