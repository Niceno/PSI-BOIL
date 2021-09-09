#include "marching_squares.h"

/******************************************************************************/
real MarchingSquares::line_density(const std::vector<LINE> & lines,
                                   const real & surf, const real & com) {
/***************************************************************************//**
*  \brief calculate total iso-line length density
*******************************************************************************/
 real length(0.);
 for(auto & l : lines)
   length += l.length();

 return length/surf;
}
