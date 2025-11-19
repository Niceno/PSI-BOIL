#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H

#include "../../Global/global_constants.h"
#include "../../Field/Vector/vector.h"
#include "../../Matter/matter.h"
#include "../../SimulationTime/simulation_time.h"
#include "lagrangian_browsing.h"
#include "particle.h"
#include <algorithm>
#include <iomanip> 

#define OFF -1 /* why is it not defined globally? */

#define TOMIYAMA_DRAG   true  /* drag coefficient. If false then Cd = 0.44 */
#define TOMIYAMA_LIFT   true  /* lift coefficient. If false then Cl = 0.15 */

/* In case of using CIT ==> SIN and VF should be set to false */
#define SIN false  /* sin profile for the interface */
#define VF  true  /* continuum volume fraction for interfacial cells */

#define CIT false  /* Corrective Interface Tracking approach. 
                     If false, then FSL will be used */

#define COALESCENCE false /*true : coalescence is considered, 
                            false: only particle rebouncing */

/////////////////
//             //
//  Lagrangian  //
//             //
/////////////////
class Lagrangian : public Scalar {

  public:
    Lagrangian(const Scalar & c,
                     Scalar * color,
                     real     lagp,  //no-use, just-for-mark
               const Scalar & cfv,  //cfvclv,  //color-function-volume-continuous-liquid-vapor
               const Vector & v,
               const Times  & t, 
               const Matter * f,
               const Matter * s);
    ~Lagrangian();

    void advance(Vector * xyz);

    void save(const char * nm, const int it);
    void load(const char * nm, const int it);

    /* add particles */
    void add(const Particle & p) {
      particles.push_back(p);
      int last_p = size() - 1; /* take the label of the last particle */
      box(last_p);
      repaint_box(last_p);
    }

    void add_cit(const Particle & p) {particles.push_back(p);}

    Particle & operator [] (int p) {return particles[p];}

    /* mathematical operators (all are inherited from Scalar, but "=" */
    const Lagrangian & operator = (const real & d)
      {for_aijk(i,j,k) val[i][j][k] = d; return *this;}

    const Domain * domain() const {return dom;}
    int size() const {return particles.size();}

    /* calculate particles's data: volume-averaged velocity, etc */
    void avg();

    /* used in particle column similation */
    void check_add_particles();

  private:

    real interface_fraction(int i, int j, int k, int p, 
                            const real & cur_val) const;

    real lagrangian_fraction(const real alfa) {
      if(lagrangian > continuous) 
        {return alfa;
        OPR(alfa);
        getchar();}
      else                      
        {return (1.0-alfa);
       	OPR(1.0-alfa);
        getchar();}          
    }                            //no-use-for-one-point-model-of-particle-tracking

          Scalar * col;
    const Scalar * cfu;  //color-function-volume-continuous-liquid-vapor
    const Vector * u;
    const Times  * time;
    const Matter * flu;  //fluid
    const Matter * sol;  //solid
    const Domain * dom; 

    void count_particles();
    void collisions();
    void couple_interface();       
    void forces();  
    void volume_averaged();  
    void modeled(Vector * xyz);  
    void merge(int p);                                   
    void merge(int pa, int pb);                                   
    void rebounce(int pa, int pb);                                   

    const real box_diam_ratio;
    const real list_diameter; /* cell size in linked_list */
    const real continuous;
    const real lagrangian;

    void box_particles() {
      for_p(p) {
        box_init(p);
      }

      /* here you fetch the velocities */
      box_velocity();

      /* correct indices for parallelization */
      for_p(p) {
        box_correct(p);
      }
    }
   
    /* repaint all particles, browsing inside particles box */
    void repaint();  

    /* repaint partilce p, browsing inside its box */
    void repaint_box(int);  

    /* filling partilce p box cells with continuous before removing p */
    void repaint_erase(int); 

    /* initialize box indices.. si, ei, etc */
    void box_init(int);

    /* fetch the velocity vector at box's corner */
    void box_velocity();

    /* correct box indices, considering domain decomposion */
    void box_correct(int);
    void box(int p) {box_init(p); box_correct(p);}

    void erase(int p) {OPR(p); OPR(size());
                       assert(p > -1);
                       assert(p < size());
                       box(p);
                       repaint_erase(p);
                       particles.erase(particles.begin() + p);
                       correct_id(p);}


    std::vector<Particle> particles;
  
    /* used for collisions */
    int * link;
    int *** cell;
    void cell_init();
    void cell_link();
    int NX_coarse, NY_coarse, NZ_coarse;
    double diam_x, diam_y, diam_z;

    Scalar p_id; /* p_id[i][j][k] gives the number of particle that
                                                is in cell(i,j,k)*/

    void correct_id(int p);

    int particle_id (int i, int j, int k) {
      int p = (int)floor(p_id[i][j][k]);
      return p;
    }
   
    /* set limit to particle velocity */
    real uvw_limit;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian.h,v 1.13 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
