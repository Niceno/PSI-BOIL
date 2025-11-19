#ifndef REALISTIC_H
#define REALISTIC_H

#include "global_constants.h"

namespace boil {
 inline bool realistic(const real val) { return std::fabs(val)<boil::zetta;  }
 inline bool realistic(const int val) { return std::abs(val)<boil::unint-10; }
}

#endif
