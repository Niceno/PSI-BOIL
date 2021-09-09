#ifndef TWOLEVELSCALAR_H
#define TWOLEVELSCALAR_H

#include "../scalar.h"

////////////////////////
//                    //
//  Two-level Scalar  //
//                    //
////////////////////////
#define for_coarsefine(l) for(int l=0; l<=1; l++)

class TwoLevelScalar {
  public:
    explicit TwoLevelScalar(const TwoLevelDomain & d,
                            const char * n = NULL) :
      coarse(d.coarse(),n),
      fine(d.fine(),n),
      levels({&coarse,&fine})
    {
    }

    Scalar & operator [] (const int l) {return *(levels[l]);}

    /* some restriction operators */
    void restrict_volume_XZ();
    void restrict_area_XZ();
    void restrict_sum_XZ();

    Scalar coarse;
    Scalar fine;

    std::array<Scalar*,2> levels;
};
#endif
