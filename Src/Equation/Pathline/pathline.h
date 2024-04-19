#ifndef PATHLINE_H
#define PATHLINE_H

#include "../../Global/global_constants.h"
#include "../../Field/Vector/vector.h"
#include "../../SimulationTime/simulation_time.h"
#include "particle.h"

////////////////
//            //
//  Pathline  //
//            //
////////////////
class Pathline {

  public:
    Pathline(const Vector & v, const Times * t);
    ~Pathline();

    void advance();

    void save(const char * nm, const int it);
    void load(const char * nm, const int it);
    void rm  (const char * nm, const int it);

    /* add particles */
    void add(const real x, const real y, const real z);

    int np() const {return npa;}

    std::vector<Particle> particles;

  protected:
    const Vector * uvw;
    const Times * time;

  private:
    void save(std::ofstream &);
    void load(std::ifstream &);
    int npa;  // number of particles 
};	

#endif
