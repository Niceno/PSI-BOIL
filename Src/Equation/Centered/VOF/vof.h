#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_realistic.h"

///////////
//       //
//  VOF  //
//       //
///////////
class VOF : public Centered {
  public:
    VOF(const Scalar & phi,
        const Scalar & f,
        const Scalar & kappa,
        const Vector & u, 
        Times & t,
        Krylov * S,
        Vector * bndclr = NULL);
    ~VOF();

    void new_time_step(){};
    void advance();
    void curvature();
    void tension(Vector * vec, const Matter matt);
    void totalvol();
    void front_minmax();
    void init(){};

    // getter for front_minmax
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};

    Vector fs;
    Vector * bndclr;
    Scalar nalpha;
    Scalar nx,ny,nz,nmag;/* normal to interface */
    Scalar adens; /* area density */
  protected:
    void advance_x();
    void advance_y();
    void advance_z();
    void bdcurv(const Scalar & g, const real & v);
    void cal_fs();
    void ext_fs();
    //void cal_fs2();
    void cal_fs3();
    void fs_bnd();
    void update_at_walls();
    void curv_HF();
    void extract_alpha();
    void insert_bc(const Scalar & g);
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_norm_cc(const Scalar & g);
    void insert_bc_norm();
    void norm_cc(const Scalar & g);
    void normalize(real & r1, real & r2, real & r3);
    real calc_v(real r1, real r2, real r3, real r4);
    real calc_alpha(real & r1, real & r2, real & r3, real & r4);
    void selectMax(const real r1, const real r2, const real r3,
                   const real r4, const real r5, const real r6,
                   const real r7, const real r8, const real r9,
                   const int i1,  const int i2,  const int i3);
    void set_iflag();
    void insert_bc_flag(ScalarInt & g, const bool b);

    real extrapolate_v(const int i, const int j, const int k,
                       const int ofx, const int ofy, const int ofz,
                       const real xp, const real yp, const real zp, real tol);

    real marching_cube_area(const int i, const int j, const int k);
#if 1
    real vel_value(const Comp m, const int i, const int j, const int k);
    real vel_correct(const int i, const int j, const int k,
                     const bool dirx, const bool diry, const bool dirz,
                     const real coef, const real mflx);
#endif
    real fext_cut(const int i, const int j, const int k, const real fval);

    void cal_adens();
    void cal_bndclr();

    void norm_cc_imin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_imax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmax(const Scalar &g, const int i,const int j, const int k);

    real alpha_val(const int i, const int j, const int k);
    real fs_val(const Comp m, const int i, const int j, const int k);
    real frontPosition(const int i, const int j, const int k, const Comp m);

    Scalar clr;     /* color function */
    Scalar kappa;        /* curvature */
    Scalar stmp,stmp2;
    Scalar fsx,fsy,fsz;
    ScalarInt iflag,iflagx,iflagy,iflagz;

    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real pi,theta;
    real dxmin,ww;
    real epsnorm;
    real phisurf;

    int nlayer, n_ext_fs;
    //int *** iflag;
};	
#endif

