#include "body.h"

/******************************************************************************/
AngleNode::AngleNode(const real an,
                     const real xp, const real yp, const real zp) {
  angle = an;
  x     = xp;
  y     = yp;
  z     = zp;
}

/******************************************************************************/
std::ostream & operator << (std::ostream &ost, const AngleNode & an) {

  ost << an.angle << " " << an.x << " " << an.y << " " << an.z;
  
  return ost;
}

/*-----------------------------------------------------------------------------+
 '$Id: anglenode.cpp,v 1.1 2008/11/29 13:20:47 niceno Exp $'/
+-----------------------------------------------------------------------------*/
