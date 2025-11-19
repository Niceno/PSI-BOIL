#ifndef REGION_H
#define REGION_H

#include <vector>
#include "../../Global/global_constants.h"

////////////////
//            //
//  Region    //
//            //
////////////////
class Region {

  public:
    Region(const int rid,
           const int cvol=0,
           const real x=0., const real y=0., const real z=0.,
           const real u=0., const real v=0., const real w=0.);

    /* get members */
    real x()      const {return m_pos[0];}
    real y()      const {return m_pos[1];}
    real z()      const {return m_pos[2];}

    real xold() const {return m_opos[0];}
    real yold() const {return m_opos[1];}
    real zold() const {return m_opos[2];}

    real u()      const {return m_vel[0];}
    real v()      const {return m_vel[1];}
    real w()      const {return m_vel[2];}

    real comu()      const {return m_comvel[0];}
    real comv()      const {return m_comvel[1];}
    real comw()      const {return m_comvel[2];}

    int cellvol() const {return m_cellvol;}
    int id() const {return m_id;}
    bool hiding() const {return m_hiding;}
    int get_tsteps_hidden() const {return m_tsteps_hiding;}

    real uvw_magn() const {return sqrt(  u()*u()
                                       + v()*v()
                                       + w()*w() );}

    int inc_hiding() {
      m_tsteps_hiding++;
      return m_tsteps_hiding;
    } 
    
    void set_oldxyz() {m_opos[0]=m_pos[0];m_opos[1]=m_pos[1];m_opos[2]=m_pos[2];}

    /* set members */
    void x(real ix) {m_pos[0]=ix;}
    void y(real iy) {m_pos[1]=iy;}
    void z(real iz) {m_pos[2]=iz;}

    void zold(real iz) {m_opos[2]=iz;}

    void xyz(real ix, real iy, real iz) {
      m_pos[0]=ix; m_pos[1]=iy; m_pos[2]=iz; }

    void oldxyz(real ox, real oy, real oz) {
      m_opos[0]=ox; m_opos[1]=oy; m_opos[2]=oz; }

    void uvw(real iu, real iv, real iw) {
      m_vel[0]=iu; m_vel[1]=iv; m_vel[2]=iw; }

    void comuvw(real iu, real iv, real iw) {
      m_comvel[0]=iu; m_comvel[1]=iv; m_comvel[2]=iw; }

    void cellvol(int icellvol) {m_cellvol = icellvol;}
    void id(int iid) {m_id=iid;}
    void set_nts_hidden(int itshid) {m_tsteps_hiding = itshid;}
    void hiding(bool ihiding) {
      m_hiding = ihiding;
      if (ihiding == false) m_tsteps_hiding = 0;
    }

  private:
    real m_pos[DIM]; 
    real m_opos[DIM];   //for position xyz at n-1
    real m_vel[DIM];
    real m_comvel[DIM];
    int m_cellvol;
    bool m_hiding;
    int m_id;
    int m_tsteps_hiding;

};

#endif
