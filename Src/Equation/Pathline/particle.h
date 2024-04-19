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
    Particle(const real r1, const real r2, const real r3);

    real x()      const {return xpos;}
    real y()      const {return ypos;}
    real z()      const {return zpos;}

    void x(real r){xpos=r;}
    void y(real r){ypos=r;}
    void z(real r){zpos=r;}

    real u()      const {return u_vel;}
    real v()      const {return v_vel;}
    real w()      const {return w_vel;}

    void u(real r){u_vel=r;}
    void v(real r){v_vel=r;}
    void w(real r){w_vel=r;}

    real vel_mag() const {return sqrt( u_vel * u_vel
                                     + v_vel * v_vel
                                     + w_vel * w_vel);}

  private:

    real xpos,ypos,zpos;
    real u_vel,v_vel,w_vel;

    /* locally numbered logical coordinates */
    int s_x, s_y, s_z;
    int e_x, e_y, e_z;

};	

#endif
