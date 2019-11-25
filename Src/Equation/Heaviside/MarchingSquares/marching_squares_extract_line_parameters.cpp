#include "marching_squares.h"

/******************************************************************************/
int MarchingSquares::extract_line_parameters(const std::vector<LINE> & lines,
                                          real & nx, real & ny, real & nalpha) {
/***************************************************************************//**
*  \brief extract the normal vector and the line constant from a line object.
*******************************************************************************/
  
  if       (lines.size()==1) {
    nx  = lines[0].normal.x;
    ny  = lines[0].normal.y;
    nalpha = lines[0].alpha;

    return 1;

  } else if(lines.size()>1) {
    real nnx(0.), nny(0.), nalp(0.), totl(0.);
    /* it can happen that there are more lines. In this case, we simplify by ave
       raging...I know it is not the best but it should do! */
    for(auto & l : lines) {
      nnx += l.normal.x*l.length(); 
      nny += l.normal.y*l.length(); 
      nalp += l.alpha*l.length(); 
      totl += l.length();
    }
    nnx /= totl;
    nny /= totl;
    nalp /= totl;

    real nsum = sqrt(nnx*nnx + nny) + boil::pico;
    nnx /= nsum;
    nny /= nsum;
    nalp /= totl;
    
    nx = nnx;
    ny = nny;
    nalpha = nalp;

    return lines.size();
  }

  /* lines is empty */
  nx = 0.;
  nx = 0.;
  nalpha = boil::unreal;

  return -1;
}
