#ifndef LevelSet_H
#define LevelSet_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

#define RCIP
#define LOCAL_OLSSON
//#define DEBUG

////////////////
//            //
//  LevelSet  //
//            //
////////////////
class LevelSet : public Centered {
  public:
    LevelSet(const Scalar & phi,
                 const Scalar & fx,
                 const Vector & u, 
                 Times & t,
                 Krylov * S);
    ~LevelSet();

    void new_time_step(){};
    void advance();
    void tension(Vector * vec, const Matter matt, Scalar & sca);
    void totalvol(const Scalar & g);
    void front_minmax();
    void update_color(Scalar & g);

    /* getter for front_minmax */
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};

    /* setter for nredist */
    void set_nredist(const int i) {
      nredist=i;
      boil::oout<<"set_nredist: nredist= "<<nredist<<boil::endl;
    }
    /* getter for nredist */
    int get_nredist() { return(nredist);};

    /* setter for cangle */
    void set_cangle(const real r) {
      cangle=r;
      boil::oout<<"set_cangle: cangle= "<<cangle<<boil::endl;
    }
    /* getter for cangle */
    real get_cangle() { return(cangle);};

  protected:
    void convection();
    void redist(const int i);
    void gradphi();
    void gradphic();
    void curv();
    void insert_bc_color(const Scalar & g);
    void insert_bc_dist(Scalar & g);
    void insert_bc_gradphi(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_kappa();
    void insert_bc_kappa2(Scalar & g);
    void insert_bc_norm();
    void bdcurv(const Scalar & g, const real & v);

    Scalar nx,ny,nz,nmag;/* normal to interface */
    Scalar kappa,stmp,dflag;
    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real pi,theta,phimin,phimax,phisurf;
    real dxmin;
    int nredist, nlayer;
    real epsnorm;
    real cangle;
};
#endif

/*-----------------------------------------------------------------------------+
 '$Id: levelset.h,v 1.2 2012/09/13 08:42:26 niceno Exp $'/
+-----------------------------------------------------------------------------*/
