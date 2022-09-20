 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <iomanip>
#  include <sstream>
#  include <stdio.h>
#  include <string>
#include "ThirdPartyHeadersEnd.h"
namespace tecplot { template <typename T> std::string ___4185( T   ___4312, int precision = 6) { std::stringstream ss; ss << std::setprecision(precision) << ___4312; return ss.str(); } } namespace tecplot { template <typename T> std::string toString(T ___4312) { return tecplot::___4185(___4312); } } template <typename T> std::string ___4185(T ___4312) { return tecplot::___4185(___4312); }
