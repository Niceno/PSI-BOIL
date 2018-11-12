#include "lookuptable.h"

/*============================================================================*/
real LookUpTable::look_up(const real val, 
                          const Column & col0, const Column & col1) const {
/*----------------------------------------------------------+
|                                                           |
|  table:   0     1     2     3     4     5     6     7     |
|           |-----|-----|-----|-----|-----|-----|-----|     |
|  seg:        1     2     3     4     5     6     7        |
|                                                           |
+----------------------------------------------------------*/

  assert(val >= table[0]             [col0]);
  assert(val <= table[table.size()-1][col0]);

  int lowseg  = 1;
  int highseg = table.size()-1;

  int currseg = (highseg+lowseg)/2;

  for(;;) {
    if(val > table[currseg][col0]) {
      lowseg=currseg;
      currseg = (int) ceil(((real)highseg+(real)lowseg)/2.0);
    } else if(val < table[currseg-1][col0]) {
      highseg=currseg;
      currseg = (int) floor(((real)highseg+(real)lowseg)/2.0);
    } else {
      const real high = table[currseg]  [col0];
      const real  low = table[currseg-1][col0];
      const real hcoef = (val -low) / (high-low);
      const real lcoef = (high-val) / (high-low);
      return   lcoef * table[currseg-1][(int)col1] 
             + hcoef * table[currseg]  [(int)col1];
    }
  }

  return -1;
}

