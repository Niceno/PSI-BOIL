#ifndef SOLIDPARTICLE_H
#define SOLIDPARTICLE_H

#include "../pathline.h"
#include "../../../Matter/matter.h"

struct ParticleProp {
    double x,y,z,diameter,density;
};


///////////////////////////////////////////////////////
//                                                   //
//  SolidParticle                                    //
//  Lagrangian particle tracking for solid           //
//  One way coupling                                 //
//  Drag and buoyancy forces are taken into account  //
//                                                   //
///////////////////////////////////////////////////////
class SolidParticle : public Pathline {
  public:
    SolidParticle(const Vector & v, const Times *t, Matter *m,
                  const real gx, const real gy, const real gz,
                  const Scalar * s1 = NULL,
                  const Scalar * s2 = NULL, const Scalar *s3 = NULL );
    ~SolidParticle() {};

    const Matter * fluid() const {return flu;}

    /* advance pathline */
    void advance();

    void save(const char * nm, const int it);
    void load(const char * nm, const int it);

    /* add particles */
    //void add_global(const real x, const real y, const real z,
    //               const real di, const real den);
    //void add_local(const real x, const real y, const real z,
    //               const real di, const real den);
    //void exchange();

    /* number of current particles */
    //int np() const {return npa;}
    //void np(const int i) {npa=i;};

    /* number of additional variables to be traced */
    //int nval() const {return nv;}

    //std::vector<Particle> particles;
    //const Scalar * s1, * s2, * s3;

  protected:
    //const Vector * uvw;
    //const Times * time;
    //const Scalar * s1, * s2, * s3;
    Matter * flu;

  private:
    void save(std::ofstream &);
    void load(std::ifstream &);
    real gravity_x, gravity_y, gravity_z;  //gravity
    //int npa;  // number of particles 
    //int nv; // number of additional variables
    //int id_serial;
    //std::vector<ParticleProp> particle_local;
};

#endif
