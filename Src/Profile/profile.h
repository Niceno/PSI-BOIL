#ifndef PROFILE_H
#define PROFILE_H

#include "../Parallel/mpi_macros.h"
#include <fstream>
#include <vector>
#include <cfloat>

#include "../Global/global_precision.h"
#include "../Global/global_split_string.h"
#include "../Parallel/Out/out.h"
#include "../Parallel/Out/print.h"

///////////////
//           //
//  Profile  //
//           //
///////////////
class Profile {
  public:
    Profile(const std::string name);
    Profile(const Profile & prof);

    real value_at(const real & x) const;

    void scale_coords(const real & s);
    void scale_values(const real & s);

    friend std::ostream & 
      operator << (std::ostream & os, const Profile & prof);

  private:
    std::vector<real> coord;  /* coordinates */
    std::vector<real> value;  /* values      */
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: profile.h,v 1.2 2014/08/06 07:45:49 sato Exp $'/
+-----------------------------------------------------------------------------*/
