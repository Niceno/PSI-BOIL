#ifndef DISTANCE_H
#define DISTANCE_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

////////////////
//            //
//  Distance  //
//            //
////////////////
class Distance : public Centered {
  public:
    Distance(const Scalar & phi, 
             const Scalar & f,   
             const Vector & u,   
             Times & t,
             Krylov * S);
    ~Distance();
	  
    void compute();

  protected:
    Matter v_fluid; /* virtual fluid for distance transport equation */

  private:
    void smooth();  /* hide it */
};	

#endif
