#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

//////////////////
//              //
//  PhaseField  //
//              //
//////////////////
class PhaseField : public Centered {
  public:
    PhaseField(const Scalar & phi, 
             const Scalar & f,
             const real & con, 
             const real & den,
             const Vector & u, 
             Times & t,
             Krylov * S);
    ~PhaseField();
	  
    void discretize(); 
    void init_smooth();
    void sharpen();
    void advance() {Centered::advance(); /*local();*/}
    void tension(Vector * vec, const Matter matt) {}

  protected:
    Scalar phi_old; /* old value */
    real   ldt;     /* local time step */

    Matter jelly;   /* virtual fluid for level set transport */
};	

#endif
