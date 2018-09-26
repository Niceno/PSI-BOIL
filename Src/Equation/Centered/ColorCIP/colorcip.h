#ifndef COLORCIP_H
#define COLORCIP_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

#define USE_TAN

////////////////
//            //
//  ColorCIP  //
//            //
////////////////
class ColorCIP : public Centered {
  public:
    ColorCIP(const Scalar & phi,
             const Scalar & f,
             const real & con, 
             const real & den,
             const Vector & u, 
             Times & t,
             Krylov * S);
    ~ColorCIP();

    void advance();
    void discretize(); 
    void tension(Vector * vec, const Matter matt);
    void totalvol();
    void front_minmax();
    void init();

    // getter for front_minmax
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};

    /* setter for ww */
    void set_ww(const real a) { 
      ww=a*dxmin;
      boil::oout<<"cipcsl2: dxmin, ww= "<<dxmin<<" "<<ww<<"\n";
    };

    /* setter for nredist */
    void set_nredist(const int i) {
      boil::oout<<"set_nredist: WARNING!!!!!!!!!!!!!!!!"<<boil::endl;
      boil::oout<<"There is no redistance function in "<<boil::endl;
      boil::oout<<"ColorCIP class."<<boil::endl;
    }

  protected:
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void curv(const Scalar & g);
    void smooth(const Scalar & g1, Scalar & g2, const int i);
    void smooth_ls(const Scalar & g1, Scalar & g2, const int i);
    void smooth_tanh(const Scalar & g1, Scalar & g2, const int i);
    void bdcurv(const Scalar & g, const real & v);
    void insert_bc(const Scalar & g);
    void insert_bc_ls(const Scalar & g);
    void insert_bc_tanh(const Scalar & g1, const Scalar & g2);


    Scalar nx,ny,nz,nmag;/* normal to interface */
    Scalar clr,clrn;     /* color function */
    Scalar gpx,gpy,gpz,gpxn,gpyn,gpzn;
    Scalar kappa,stmp,ssp,diag,rsdl;
    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real pi,tanfac,theta;
    real dxmin,ww;
    real epsnorm;

    int *** iflag;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: colorcip.h,v 1.6 2009/11/12 12:24:07 sato Exp $'/
+-----------------------------------------------------------------------------*/

