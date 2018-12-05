#ifndef REALISTIC_H
#define REALISTIC_H

#include "global_constants.h"

namespace boil {
 inline bool realistic(const real val) { return val<boil::zetta; }
}

#endif
