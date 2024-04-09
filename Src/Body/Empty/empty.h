#ifndef EMPTY_H
#define EMPTY_H

#include "../body.h"

/////////////
//         //
//  EMPTY  //
//         //
/////////////
class Empty : public Body {
  public:
    Empty() {}

    bool cut(int i, int j, int k) const {return false;}
    bool off(int i, int j, int k) const {return false;}
    bool on (int i, int j, int k) const {return true;}

    bool cut(const Comp & m, int i, int j, int k) const {return false;}
    bool off(const Comp & m, int i, int j, int k) const {return false;}
    bool on (const Comp & m, int i, int j, int k) const {return true;}

    real fSw(const Comp & m, int i, int j, int k) const {return 1.0;}
    real fSe(const Comp & m, int i, int j, int k) const {return 1.0;}
    real fSs(const Comp & m, int i, int j, int k) const {return 1.0;}
    real fSn(const Comp & m, int i, int j, int k) const {return 1.0;}
    real fSb(const Comp & m, int i, int j, int k) const {return 1.0;}
    real fSt(const Comp & m, int i, int j, int k) const {return 1.0;}

};

#endif
