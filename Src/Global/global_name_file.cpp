#include "global_name_file.h"

/******************************************************************************/
const std::string name_file(const char * gname, const char * ext, 
                            const int ts,    /* time step or grid level */
                            const int cpu) { /* (fake) processor i.d. */

  /* start with the given name */
  std::ostringstream numb;

  numb << gname;

  /* add processor number */
  if( cpu > -1 ) {
    numb << "_p";
    numb.fill('0');
    if( cpu < 1000) {
      numb.width(3);
    } else {
      numb.width(4);
    }
    numb << cpu;
  }

  /* add time step (or level) */
  if( ts > -1 ) {
    numb << "_";
    numb.fill('0');
    numb.width(6);
    numb << ts;
  }

  /* add extension */
  numb << ext;

  /* finalize the file name */
  return numb.str();
}
