#include <iostream>
#include <sstream>

#ifndef NAME_FILE_H
#define NAME_FILE_H

const std::string name_file(const char * gname, const char * ext, 
                            const int ts,       /* time step or grid level */
                            const int cpu=-1);  /* (fake) processor i.d. */

#endif

/*-----------------------------------------------------------------------------+
 '$Id: global_name_file.h,v 1.4 2008/11/17 19:23:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
