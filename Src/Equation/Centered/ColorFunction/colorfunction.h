#ifndef COLORFUNCTION_H
#define COLORFUNCTION_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

/////////////////////
//                 //
//  ColorFunction  //
//                 //
/////////////////////
class ColorFunction : public Centered {
  public:
    ColorFunction(const Scalar & phi, 
             const Scalar & f,
             const real & con, 
             const real & den,
             const Vector & u, 
             Times & t,
             Krylov * S);
    ~ColorFunction();
	  
    void discretize(); 
    void init_smooth();
    void sharpen();
    void advance() {Centered::advance(); /*local();*/}
    void tension(Vector * vec, const Matter matt);
    void totalvol();
    void front_minmax();
    void front_minmax( Range<real> xr
                     , Range<real> yr
                     , Range<real> zr );
    /* getter for front_minmax */
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};
    /* getter for clrsum1 & clrsum2 */
    real get_clrsum1() { return (clrsum1); };
    real get_clrsum2() { return (clrsum2); };

  protected:
    Scalar nx,ny,nz;/* normal to interface */
    Scalar ndiv;
    Scalar surf;
    Scalar phi_old; /* old value */
    real   ldt;     /* local time step */
    real   xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real   clrsum1,clrsum2;
    Matter v_fluid; /* virtual fluid for level set transport */
};	

#endif
