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

    //! Is empty
    virtual bool is_empty() { return true; }

    virtual bool cut(int i, int j, int k) const {return false;}
    virtual bool off(int i, int j, int k) const {return false;}
    virtual bool on (int i, int j, int k) const {return true;}

    virtual bool cut(const Comp & m, int i, int j, int k) const {return false;}
    virtual bool off(const Comp & m, int i, int j, int k) const {return false;}
    virtual bool on (const Comp & m, int i, int j, int k) const {return true;}

    /* only for pressure equation */
    virtual bool cut_p(int i, int j, int k) const { return false; };
    virtual bool off_p(int i, int j, int k) const { return false; };
    virtual bool on_p (int i, int j, int k) const { return true; };

    virtual real fSw(const Comp & m, int i, int j, int k) const {return 1.0;}
    virtual real fSe(const Comp & m, int i, int j, int k) const {return 1.0;}
    virtual real fSs(const Comp & m, int i, int j, int k) const {return 1.0;}
    virtual real fSn(const Comp & m, int i, int j, int k) const {return 1.0;}
    virtual real fSb(const Comp & m, int i, int j, int k) const {return 1.0;}
    virtual real fSt(const Comp & m, int i, int j, int k) const {return 1.0;}

    virtual real fdxw(int i, int j, int k) const { return 1.0; } 
    virtual real fdxe(int i, int j, int k) const { return 1.0; } 
    virtual real fdys(int i, int j, int k) const { return 1.0; } 
    virtual real fdyn(int i, int j, int k) const { return 1.0; } 
    virtual real fdzb(int i, int j, int k) const { return 1.0; } 
    virtual real fdzt(int i, int j, int k) const { return 1.0; } 

};

#endif
