#include "rect.h"

Rect
operator||( const Rect & rectA, const Rect & rectB )
{
  Rect rect;
  rect.top = ( rectA.top > rectB.top ) ? rectA.top : rectB.top;
  rect.left = (rectA.left < rectB.left) ? rectA.left : rectB.left;
  if ( rectA.left + rectA.width > rectB.left + rectB.width )
    rect.width = rectA.left + rectA.width - rect.left;
  else
    rect.width = rectB.left + rectB.width - rect.left;
  if ( rectA.top - rectA.height < rectB.top - rectB.height )
    rect.height = rect.top - ( rectA.top - rectA.height );
  else
    rect.height = rect.top - ( rectB.top - rectB.height );
  return rect;
}

/*-----------------------------------------------------------------------------+
 '$Id: rect.cpp,v 1.2 2008/10/21 11:54:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
