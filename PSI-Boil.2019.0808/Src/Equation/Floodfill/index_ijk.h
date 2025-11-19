#ifndef INDEX_IJK_H
#define INDEX_IJK_H
#include "../../Field/Scalar/scalar.h"

//////////////
//          //
//  Indijk  //
//          //
//////////////
class Index_ijk {
  private:
    int xp[3];

  public:
    Index_ijk()                       {xp[0] = 0;  xp[1] = 0;  xp[2] = 0;}
    Index_ijk(int xa, int ya, int za) {xp[0] = xa; xp[1] = ya; xp[2] = za;}

    void set_to(int xa, int ya, int za) {xp[0]=xa;   xp[1]=ya;   xp[2]=za;}
    void set_to(const Index_ijk & o) {for(int i=0; i<3; i++) xp[i] = o.xp[i];}

    int i() const {return xp[0];}
    int j() const {return xp[1];}
    int k() const {return xp[2];}
};

#endif
