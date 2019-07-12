#ifndef VOF_H
#define VOF_H

#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_realistic.h"
#include "../../Heaviside/heaviside.h"
#include "../../Topology/topology.h"

#define IB

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
    void advance(Scalar & sca);
    void curvature();
    void ancillary(); /* calcs ancillary params such as adens w/o advance */
    void tension(Vector * vec, const Matter matt);
    void tension(Vector * vec, const Matter matt, const Scalar & scp);
    void totalvol();
    void front_minmax();
    void init(){};

    void cal_liq_vel(Vector * umass, Vector * uliq);
    void ext_vel(Scalar & sca, const Scalar & eflag, const int sgn);
    void smooth(const Scalar & sca, Scalar & scb, const int itnum);

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

    /* setter for cangle */
    void set_cangle(const real r) {
      cangle=r/180.0*acos(-1.0);
      boil::oout<<"set_cangle: cangle= "<<r<<"\n";
      boil::oout<<"#############################################\n";
      boil::oout<<"# WARNING!!! cangle is not implemented yet. #\n";
      boil::oout<<"#############################################\n";
    }
    /* getter for cangle */
    real get_cangle() { return(cangle/acos(-1.0)*180.0);};

    /* setter for cangle */
    void set_limit_color(const bool b) {
      limit_color=b;
      boil::oout<<"set_limit_color= "<<b<<"\n";
    }
    /* getter for cangle */
    bool get_limit_color() { return(limit_color);};


    Vector fs;
    Vector * bndclr;
    Topology  topo;

    Scalar nalpha;
    Scalar nx,ny,nz;/* normal to interface */
    Scalar adens;
    Scalar mx,my,mz;/* normal to interface, in real space */

  protected:
    void advance_x(Scalar & sca);
    void advance_y(Scalar & sca);
    void advance_z(Scalar & sca);
    void bdcurv();
    void cal_fs3();
    void cal_fs_interp();
    void curv_HF();
    void curv_smooth();
    void extract_alpha();
    void fs_bnd();
    void standardized_norm_vect();
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void ib_norm(const Scalar & g);
    void ib_norm_cal(const int cc, const int i, const int j, const int k);
    void insert_bc(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_gradphi(const Scalar & g);
    void insert_bc_norm_cc(const Scalar & g);
    void insert_bc_norm();
    void nib(const real & n1, const real & n2, const real & n3
           , const real & n4, const real & n5, const real & n6
           , real r[]);
    void norm_cc(const Scalar & g);
    void norm_young(const Scalar & g);
    void normalize(real & r1, real & r2, real & r3);
    void nwall(const Scalar & g
             , const real & r1, const real & r2, const real & r3
             , const int & i1, const int & i2, const int &i3
             , real r[] );
    void true_norm_vect();
    void update_at_walls();
    void wall_norm(const Scalar & sca);
    real kappa_ave(const real r1, const real r2);
    real kappa_ave(const real r1, const real r2, const int i1, const int i2);

    real calc_v(real r1, real r2, real r3, real r4);
    real calc_alpha(const real r1, const real r2, const real r3, const real r4);
    real calc_flux(const real g, real c, const real nx, const real ny, const real nz);

    void selectMax(const real r1, const real r2, const real r3,
                   const real r4, const real r5, const real r6,
                   const real r7, const real r8, const real r9,
                   const int i1,  const int i2,  const int i3);
    void set_iflag();
    void insert_bc_flag(ScalarInt & g, const bool b);

    void cal_adens();
    void cal_bndclr();

    void set_adens(const Scalar & newadens) {
      for_aijk(i,j,k)
        adens[i][j][k] = newadens[i][j][k];
    }

#if 1
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


    real extrapolate_v(const int i, const int j, const int k,
                       const int ofx, const int ofy, const int ofz,
                       const real xp, const real yp, const real zp);

    void norm_cc_imin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_imax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmax(const Scalar &g, const int i,const int j, const int k);
    void vf_limiter();   

    /* elvira functions */
    void norm_elvira(Scalar & sca);
    void norm_elvira(int i,int j,int k,
                     real valcc,real valmc,real valpc,real valcm,
                     real valcp,real valmm,real valpm,real valmp,real valpp);
    void normalize_elvira(const real m, const real sig, 
                          real & nnx, real & nny, real & nnz);
    real elvira_l2(const real alp,const real nnx,const real nny,const real nnz,
                   const real valcc,const real valmc,const real valpc,
                   const real valcm,const real valcp,const real valmm,
                   const real valpm,const real valmp,const real valpp);
    real ext_v(const real xp, const real yp, const real zp, 
               const real vv1, const real vv2, const real vv3,
               const real vn1, const real vn2, const real vn3,
               const real denom, const real alp);
    void norm_cc(Scalar & sca,int i,int j,int k);

    real alpha_val(const real c, const real nnx, const real nny, const real nnz);
    real fs_val(const Comp m, const int i, const int j, const int k);
    real frontPosition(const int i, const int j, const int k, const Comp m);

    Scalar clr;     /* color function */
    Scalar kappa;        /* curvature */
    Scalar stmp, stmp2, stmp3;
    ScalarInt iflag,iflagx;

    Scalar utx, uty, utz, unliq;

    real rhol, rhov; /* densities for velocity correction */
    const Matter * mixt() const {return mixture;}
    Matter * mixture;

    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real theta;
    real epsnorm;
    real phisurf;
    real tol_wall, tol_flux, tol_ext, flux_cfl;
    real ww, dxmin;
    bool iminp, imaxp, jminp, jmaxp, kminp, kmaxp; // true = periodic
    bool iminw, imaxw, jminw, jmaxw, kminw, kmaxw; // true = wall
    bool iminc, imaxc, jminc, jmaxc, kminc, kmaxc; // true = cut-stencil
    bool ifull, jfull, kfull; // true = not a dummy direction
    bool limit_color;

    Heaviside heavi;

    int nlayer;
    int curv_method;
    real cangle;
};	
#endif

