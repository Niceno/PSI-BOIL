#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "../../Global/global_constants.h"
#include "../../Ravioli/comp.h"

////////////////
//            //
//  Particle  //
//            //
////////////////
class Particle {

  public:
    Particle(const real r1, const real r2, const real r3, const int nval);

    real x() const {return xpos;}
    real y() const {return ypos;}
    real z() const {return zpos;}

    void x(real r){xpos=r;}
    void y(real r){ypos=r;}
    void z(real r){zpos=r;}

    real u() const {return u_vel;}
    real v() const {return v_vel;}
    real w() const {return w_vel;}

    void u(real r){u_vel=r;}
    void v(real r){v_vel=r;}
    void w(real r){w_vel=r;}

    real vel_mag() const {return sqrt( u_vel * u_vel
                                     + v_vel * v_vel
                                     + w_vel * w_vel);}

    /* additional scalar to be traced */

    real sval(const int i) const {return sca[i];}
    void sval(const int i, const real r){sca[i]=r;}

    int id() const {return iid;}
    void id(int i){iid=i;}

  private:

    real xpos,ypos,zpos;
    real u_vel,v_vel,w_vel;
    int  iid;
    std::vector<double> sca;

};

#endif
