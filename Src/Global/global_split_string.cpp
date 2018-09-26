#include "global_split_string.h"

namespace boil {

/******************************************************************************/
const std::vector<std::string> split_string(std::string input, 
                                            const char * delim) {


  std::vector<std::string> result;

  int cut_at;

  while( (cut_at = input.find_first_of(delim)) != input.npos )
   {
    if(cut_at > 0)
     {
      result.push_back(input.substr(0,cut_at));
     }
    input = input.substr(cut_at+1);
   }

  if(input.length() > 0)
   {
    result.push_back(input);
   }

  return result;
}

} /* boil */

/*-----------------------------------------------------------------------------+
 '$Id: global_split_string.cpp,v 1.1 2009/05/07 17:56:44 niceno Exp $'/
+-----------------------------------------------------------------------------*/
