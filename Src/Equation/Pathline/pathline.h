#ifndef PATHLINE_H
#define PATHLINE_H

#include "../../Global/global_constants.h"
#include "../../Field/Vector/vector.h"
#include "../../SimulationTime/simulation_time.h"
#include "particle.h"

struct Coordinate {
    double x,y,z;
};

////////////////
//            //
//  Pathline  //
//            //
////////////////
class Pathline {

  public:
    Pathline(const Vector & v, const Times *t, const Scalar * s1 = NULL,
             const Scalar * s2 = NULL, const Scalar *s3 = NULL );
    ~Pathline();

    /* advance pathline */
    void advance();
    /* initialize */
    void init();

    void save(const char * nm, const int it);
    void load(const char * nm, const int it);
    void rm  (const char * nm, const int it);

    /* add particles */
    void add_global(const real x, const real y, const real z);
    void add_local(const real x, const real y, const real z);
    void exchange();

    /* number of current particles */
    int np() const {return npa;}
    void np(const int i) {npa=i;};

    /* number of additional variables to be traced */
    int nval() const {return nv;}

    std::vector<Particle> particles;
    const Scalar * s1, * s2, * s3;

  protected:
    const Vector * uvw;
    const Times * time;
    //const Scalar * s1, * s2, * s3;

  private:
    void save(std::ofstream &);
    void load(std::ifstream &);
    int npa;  // number of particles 
    int nv; // number of additional variables
    int id_serial;
    std::vector<Coordinate> particle_local;
};	

#endif
