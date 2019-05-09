#ifndef VOF_H
#define VOF_H

#include <cmath>
#include <list>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_realistic.h"
#include "../../Heaviside/heaviside.h"
#include "../../Topology/topology.h"

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
        Vector * bndclr = NULL,
        Matter * flu = NULL);
    ~VOF();

    void new_time_step(){};
    void advance();
    void ancillary(); /* calcs ancillary params such as adens w/o advance */
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

    // getter and setter for wall value tolerance
    real get_tol_wall() { return tol_wall; }
    void set_tol_wall(real tolnew) {
      tol_wall = tolnew;
      boil::oout<<"VOF: New wall value tolerance: "<<tol_wall<<boil::endl;
      return;
    }

    // getter and setter for curv_method
    int get_curv_method() {return curv_method;}
    void set_curv_method(int i) {
      curv_method=i;
      if(i==0){
        boil::oout<<"VOF: height function is used for curvature calculation.\n";
      } else if(i==1){
        boil::oout<<"VOF: smoothed VOF is used for curvature calculation.\n";
      } else {
        boil::oout<<"method should be 0 or 1.\n";
        boil::oout<<"0 for height function.\n";
        boil::oout<<"1 for smoothed VOF.\n";
        exit(0);
      }
    }

    // getter and setter for flux iteration tolerance
    real get_tol_flux() { return tol_flux; }
    void set_tol_flux(real tolnew) { 
      tol_flux = tolnew;
      boil::oout<<"VOF: New flux iteration tolerance: "<<tol_flux<<boil::endl;
      return;
    }

    // getter and setter for flux iteration number
    real get_iter_flux() { return maxiter; }
    void set_iter_flux(int iternew) {
      maxiter = iternew;
      boil::oout<<"VOF: New flux iteration number: "<<maxiter<<boil::endl;
      return;
    }

    // getter/setter for flux cfl
    real get_flux_cfl() { return flux_cfl; }
    void set_flux_cfl(real cflnew) {
      flux_cfl = cflnew;
      boil::oout<<"VOF: New flux CFL number: "<<flux_cfl<<boil::endl;
      return;
    }


    void cal_liq_vel();
    void set_adens(const Scalar & newadens) {
      for_aijk(i,j,k)
        adens[i][j][k] = newadens[i][j][k];
    }
    void sharpen();

    Scalar unliq; /* normal component of liquid velocity */
    Scalar utliq, utx, uty, utz; /* tangential component of liquid velocity */
    Vector uliq; /* liquid velocity */
    Vector * bndclr;

    Topology topo;
  protected:
    void advance_x();
    void advance_y();
    void advance_z();
    void bdcurv(const Scalar & g, const real & v);
    void cal_fs3();
    void cal_fs_interp();
    void ext_vel(Scalar & sca, const Scalar & eflag, const int sgn);
    void fs_bnd();
    void update_at_walls();
    void curv_HF();
    void curv_HF_ext();
    void curv_smooth();
    real kappa_ave(const real r1, const real r2);
    void smooth(const Scalar & sca, Scalar & scb, const int itnum);
    void extract_alpha();
    void true_norm_vect();
    void insert_bc(const Scalar & g);
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_norm_cc(const Scalar & g);
    void insert_bc_norm();
    void norm_cc(const Scalar & g);
    void norm_young(const Scalar & g);
    void normalize(real & r1, real & r2, real & r3);
    real calc_v(real r1, real r2, real r3, real r4);
    real calc_alpha(const real r1, const real r2, const real r3, const real r4);
    real calc_flux(const real g, real c, const real nx, const real ny, const real nz);

    real calc_diabatic_flux(const real jv, const real gliq, const real ggas,
                            real phiup,
                            const real nx, const real ny, const real nz);
    real calc_diabatic_flux(real & jv, const real gliq, const real ggas,
                            const real dxrat, real phiup, real phidn,
                            const real nxup, const real nyup, const real nzup,
                            const real nxdn, const real nydn, const real nzdn);
/* underdevelopment */
#if 0
    real calc_diabatic_flux(const real jv, const real gliq, 
                            const real ggasup,const real ggasdn,
                            const real dxrat, real phiup, real phidn,
                            const real nxup, const real nyup, const real nzup,
                            const real nxdn, const real nydn, const real nzdn);
#endif

    real iterate_flux(const real jv, const real dirsgn, const real cflrat,
                      const real alphaliq, const real alphagas,
                      const real vm1, const real vm2, const real vm3);
    real iterate_flux(const real jv, const real dirsgnup, const real dirsgndn,
                      const real cflrat, const real alphaup,
                      const real vm1up, const real vm2up, const real vm3up,
                      const real alphadn,
                      const real vm1dn, const real vm2dn, const real vm3dn,
                      const real x0start, const real x2start);

    real calc_v_iter(const real alpha, const real vm1, const real vm2,
                     const real vm3, const real absg, const real dirsgn);
 
    void selectMax(const real r1, const real r2, const real r3,
                   const real r4, const real r5, const real r6,
                   const real r7, const real r8, const real r9,
                   const int i1,  const int i2,  const int i3);
    void set_iflag();
    void superpose();
    void insert_bc_flag(ScalarInt & g, const bool b);

    real extrapolate_v(const int i, const int j, const int k,
                       const int ofx, const int ofy, const int ofz,
                       const real xp, const real yp, const real zp);

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
    void vf_limiter();   

    real alpha_val(const int i, const int j, const int k);
    real fs_val(const Comp m, const int i, const int j, const int k);
    real frontPosition(const int i, const int j, const int k, const Comp m);

    Scalar clr;     /* color function */

#if 0
    /* adensgeom stuff */
    typedef struct {
      real x,y,z;
    } XYZ;
    real calc_area(const std::vector<XYZ> &vect, const XYZ norm);
    XYZ PlusXYZ(const XYZ p1, const XYZ p2);
    real DotProduct(const XYZ p1, const XYZ p2);
    XYZ CrossProduct(const XYZ p1, const XYZ p2);
    void cal_adens_geom(Scalar & eval);
    Scalar adensgeom; /* area density (geometric) */

    /* mc stuff */
    typedef struct {
       XYZ p[8];
       real val[8];
    } GRIDCELL;
    XYZ VertexInterpVOF(real isolevel, XYZ p1, XYZ p2, real valp1, real valp2);
    real PolygoniseVOF(GRIDCELL grid, real isolevel);
    typedef struct {
       XYZ p[3];
       real area(){
         real x1 = p[1].x-p[0].x;
         real y1 = p[1].y-p[0].y;
         real z1 = p[1].z-p[0].z;
         real x2 = p[2].x-p[0].x;
         real y2 = p[2].y-p[0].y;
         real z2 = p[2].z-p[0].z;
         real area=  (y1*z2-z1*y2)*(y1*z2-z1*y2)
                    +(z1*x2-x1*z2)*(z1*x2-x1*z2)
                    +(x1*y2-y1*x2)*(x1*y2-y1*x2);
         area = 0.5 * sqrt(area);
         return(area);
       }
    } TRIANGLE;
#endif

    Scalar nalpha;
    Vector fs;
    Scalar nx,ny,nz,nmag;/* normal to interface */
    Scalar mx,my,mz;/* normal to interface, in real space */
    Scalar adens; /* area density */

    Heaviside heavi;

    Scalar kappa;        /* curvature */
    Scalar stmp,stmp2,stmp3;
    ScalarInt iflag,iflagx,iflagy,iflagz;
    Vector sosflux,fluxmax;
    Scalar stmp4,stmp5,stmp6;

    real rhol, rhov; /* densities for velocity correction */
    const Matter * mixt() const {return mixture;}
    Matter * mixture;

    real tol_wall, tol_flux, tol_ext, flux_cfl;
    int maxiter;

    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real pi,theta;
    real dxmin,ww;
    real epsnorm;
    real kappa_non_cal;
    bool iminp, imaxp, jminp, jmaxp, kminp, kmaxp; // periodic = true

    real phisurf;
#if 0
    real f_w, f_e, f_t, f_b, f_n, f_s;
#endif

    int nlayer, n_ext_fs;
    int curv_method;
};	
#endif

