#include "ms_axisym.h"

/******************************************************************************/
real MSaxisym::line_density(const std::vector<LINE> & lines,
                            const real & surf, const real & com) {
/***************************************************************************//**
*  \brief calculate total iso-line length density
*******************************************************************************/
 real length(0.);
 for(auto & l : lines)
   length += l.length()*l.CoM().x;

 return length/(surf*com);
}
