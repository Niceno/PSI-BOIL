#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "../../Global/global_constants.h"
#include "../../Ravioli/comp.h"
#include "../../Ravioli/diameter.h"
#include "position.h"

////////////////
//            //
//  Particle  //
//            //
////////////////
class Particle {

  public:
    Particle(const Position & pos, 
             const Diameter & dia);
    Particle(const Position & pos, 
             const Diameter & dia, 
             const Position & vel);

    real x()      const {return position.x();}
    real y()      const {return position.y();}
    real z()      const {return position.z();}

    real u()      const {return velocity.x();}
    real v()      const {return velocity.y();}
    real w()      const {return velocity.z();}

    real d()      const {return diameter.d();}
    real area()   const {return are;}
    real volume() const {return vol;}

    real & xyz     (const Comp & m) {return position.xyz(m);}
    real & uvw     (const Comp & m) {return velocity.xyz(m);}
    real   uvw_magn() const {return sqrt(  velocity.x() * velocity.x()
                                         + velocity.y() * velocity.y()
                                         + velocity.z() * velocity.z());}

    /* set locally numbered logical coordinates */
    void si(int sx) {s_x = sx;}
    void sj(int sy) {s_y = sy;}
    void sk(int sz) {s_z = sz;}
    void ei(int ex) {e_x = ex;}
    void ej(int ey) {e_y = ey;}
    void ek(int ez) {e_z = ez;}
    /* access locally numbered logical coordinates */
    int si() const {return s_x;}
    int sj() const {return s_y;}
    int sk() const {return s_z;}
    int ei() const {return e_x;}
    int ej() const {return e_y;}
    int ek() const {return e_z;}

    real box_dx() const {return boxdim[~Comp::i()];}
    real box_dy() const {return boxdim[~Comp::j()];}
    real box_dz() const {return boxdim[~Comp::k()];}
    void box_dx(real v) {boxdim[~Comp::i()] = v;}
    void box_dy(real v) {boxdim[~Comp::j()] = v;}
    void box_dz(real v) {boxdim[~Comp::k()] = v;}

    real box_u(int n) const {assert(n < 8); return boxuvw[n][~Comp::i()];}
    real box_v(int n) const {assert(n < 8); return boxuvw[n][~Comp::j()];}
    real box_w(int n) const {assert(n < 8); return boxuvw[n][~Comp::k()];}
    real box_uvw(int n, const Comp m) const {
         assert(n < 8);return boxuvw[n][~m];}

    void box_u(int n, real v) {assert(n < 8); boxuvw[n][~Comp::i()] = v;}
    void box_v(int n, real v) {assert(n < 8); boxuvw[n][~Comp::j()] = v;}
    void box_w(int n, real v) {assert(n < 8); boxuvw[n][~Comp::k()] = v;}
    void box_uvw(int n, const Comp m, real v) {assert(n < 8); boxuvw[n][~m] = v;}

    real vol_rho(const Comp m) const {return volrho[~m];}
    void vol_rho(const Comp m, real v) {volrho[~m] = v;}

    real vol_u() const {return voluvw[~Comp::i()];}
    real vol_v() const {return voluvw[~Comp::j()];}
    real vol_w() const {return voluvw[~Comp::k()];}
    real vol_uvw(const Comp m) const {return voluvw[~m];}

    void vol_u(real v) {voluvw[~Comp::i()] = v;}
    void vol_v(real v) {voluvw[~Comp::j()] = v;}
    void vol_w(real v) {voluvw[~Comp::k()] = v;}
    void vol_uvw(const Comp m, real v) {voluvw[~m] = v;}

    real uvwc_old(const Comp m) const {return u_old[~m];}
    void uvwc_old(const Comp m, real v) {u_old[~m] = v;}

    void gradient(int dir,const Comp m, real v) {grad[dir][~m] = v;}
    real gradient(int dir,const Comp m) const {return grad[dir][~m];}

  private:
    Position position;
    Position velocity;
    Diameter diameter;
    real     are;
    real     vol;
    real     mas;

    /* locally numbered logical coordinates */
    int s_x, s_y, s_z;
    int e_x, e_y, e_z;

    /* boxed values */
    std::vector< real >              boxdim;
    std::vector< std::vector<real> > boxuvw;

    /* volume averaged value */
    std::vector<real> volrho; 
    std::vector<real> voluvw; 

    std::vector<real>  u_old; 
    std::vector<real>  uvw_c; 

    std::vector< std::vector<real> > grad;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: particle.h,v 1.7 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
