#ifndef PRINT_H
#define PRINT_H

#include "out.h"    

/* All (processors) Print */
#define APR(x) boil::aout << "PROC: "   << boil::cart.iam() \
                          << ", FILE: " << __FILE__         \
                          << ", LINE: " << __LINE__         \
                          << ", " << #x " = " << x          \
                          << '\n';         

#define AMS(x) boil::aout << "PROC: "   << boil::cart.iam() \
                          << ", FILE: " << __FILE__         \
                          << ", LINE: " << __LINE__         \
                          << ", " << #x                     \
                          << '\n';         

/* One (processor) Prints */
#define OPR(x) boil::oout << "FILE: "   << __FILE__ \
                          << ", LINE: " << __LINE__ \
                          << ", " << #x " = " << x  \
                          << '\n';         

#define OMS(x) boil::oout << "FILE: "   << __FILE__ \
                          << ", LINE: " << __LINE__ \
                          << ", " << #x             \
                          << '\n';         

#endif
