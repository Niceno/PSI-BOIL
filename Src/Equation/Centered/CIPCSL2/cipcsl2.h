#ifndef CIPCSL2_H
#define CIPCSL2_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Global/global_constants.h"
#include "../../../Global/global_realistic.h"
#include "../../Topology/topology.h"
#include "../../Heaviside/MarchingCubes/marching_cubes.h"

#define RCIP
//#define LOCAL_OLSSON
#define IB
//#define DEBUG

//////////////
//          //
//  Scheme  //
//          //
//////////////
class Scheme {
  public:
    Scheme(const Scalar & phi) {
      sigx = phi.shape();
      sigx.allocate(phi.ni(),   phi.nj()+1, phi.nk()+1); 
                                sigx.oy(1); sigx.oz(1);

      sigy = phi.shape();
      sigy.allocate(phi.ni()+1, phi.nj(),   phi.nk()+1);
                    sigy.ox(1);             sigy.oz(1);

      sigz = phi.shape(); 
      sigz.allocate(phi.ni()+1, phi.nj()+1, phi.nk()  );
                    sigz.ox(1); sigz.oy(1);

      f    = phi.shape();
      f   .allocate(phi.ni()+1, phi.nj()+1, phi.nk()+1);
                    f.ox(1);    f.oy(1);    f.oz(1);
    }

    void bdcond_i(const Scalar & sca);
    void bdcond_j(const Scalar & sca);
    void bdcond_k(const Scalar & sca);
    void bdcond_f(const Scalar & sca);

    Scalar sigx; /* "scalar " */
    Scalar sigy; /* "scalar " */
    Scalar sigz; /* "scalar " */
    Scalar f;    /* nodal values */

  protected:
    Scheme() {}
};

///////////////
//           //
//  CIPCSL2  //
//           //
///////////////
class CIPCSL2 : public Centered {
  public:
    CIPCSL2(const Scalar & phi,
            const Scalar & fx,
            const Scalar & kappa,
            const Vector & u, 
            Times & t,
            Krylov * S);
    ~CIPCSL2();

    void new_time_step(){};
    void advance();
    void bdcond(const Scalar & sca);
    void curvature();
    void ib_bdcond(const Scalar & sca);
    void ib_ext_scalar(Scalar & g);
    void tension(Vector * vec, const Matter matt);
    void tension(Vector * vec, const Matter matt, Scalar & sca);
    void totalvol();
    void totalvol( Range<real> xr
                 , Range<real> yr
                 , Range<real> zr );
    void front_minmax();
    void front_minmax( Range<real> xr
                     , Range<real> yr
                     , Range<real> zr );
    void init();
    void update_node(Scalar & g);
    void save(char *, const int);
    void load(char *, const int);
    void rm(char *, const int);
    void range(){
      boil::oout<<"Range of color function: "<<clr.min()<<" "<<clr.max()<<"\n";
    }
    void bnd_extract( const Dir d, real *** cp ) const;
    void bnd_insert ( const Dir d, real **  cp );

    /* getter for front_minmax */
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
      nredist=i;
      boil::oout<<"set_nredist: nredist= "<<nredist<<boil::endl;
    }
    /* getter for nredist */
    int get_nredist() { return(nredist);};

    /* setter for itsharpen */
    void set_itsharpen(const int i) {
      itsharpen=i;
      boil::oout<<"set_itsharpen: itsharpen= "<<itsharpen<<boil::endl;
    }
    /* getter for itsharpen */  
    int get_itsharpen() { return(itsharpen);};

    /* setter for localSharpen */
    void set_localSharpen() { localSharpen=true;};
    void set_globalSharpen() { localSharpen=false;};
    /* getter for localSharpen */
    bool get_localSharpen() {return localSharpen;};

    /* setter for use_dist_for_kappa */
    void set_use_dist_for_kappa(bool b) { use_dist_for_kappa=b;
      boil::oout<<"use_dist_for_kappa= "<<use_dist_for_kappa<<"\n";};
    /* getter for use_dist_for_kappa */
    bool get_use_dist_for_kappa() {return use_dist_for_kappa;};

    /* setter for cangle */
    void set_cangle(const real r) {
      cangle=r/180.0*acos(-1.0);
      boil::oout<<"set_cangle: cangle= "<<r<<"\n";
    }
    /* getter for cangle */
    real get_cangle() { return(cangle/acos(-1.0)*180.0);};

    /* setter for eps_st */
    void set_eps_st(const real r) {
      eps_st=r;
      boil::oout<<"set_eps_st: eps_st= "<<eps_st<<boil::endl;
    }
    /* getter for eps_st */
    real get_eps_st() { return(eps_st);};

    /* setter for itsmear */
    void set_itsmear(const int i) {
      itsmear=i;
      boil::oout<<"itsmear: itsmear= "<<itsmear<<boil::endl;
    }
    /* getter for itsmear */
    real get_itsmear() { return(itsmear);};

    /* getter for sum_outlet */
    real get_sum_outlet() { return (sum_outlet); };

    /* getter for sum_outletm */
    real get_sum_outletm() { return (sum_outletm); };

    /* getter for clrsum1 */
    real get_clrsum1() { return (clrsum1); };

    /* getter for clrsum2 */
    real get_clrsum2() { return (clrsum2); };

    /* min and max of color function in fluid domain */
    real minval() {return minclr;}
    real maxval() {return maxclr;}
    void color_minmax(); 

    Topology topo;
    void ancillary();

  protected:
    void bnd_wall_kappa();
    void convection();
    void curv(const Scalar & g);
    void curv_interface();
    void curv_interface_ext();
    void bdcurv_interface();
    void bdcurv_interface_ext();
    void distfunc(const Scalar & g, const int i);
    void ext_sca(Scalar & g);
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void set_iflag();
    void set_wflag();
    void set_wflag2();
    void normalize(real & r1, real & r2, real & r3);
    void smear(Scalar & g);
    bool sharpen(Scalar & g, const real e, const int i, const bool b);
    void insert_bc_alp(const Scalar & g);
    void insert_bc_dist(Scalar & g);
    void insert_bc_flag(ScalarInt & g, const bool b);
    void insert_bc_gradphi(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_kappa(Scalar & g);
    void insert_bc_norm();
    void redist(const bool b);
    void set_alp();
    void debugout(const char *);
    void bdcurv(const Scalar & g);
    void ib_norm(const Scalar & g);
    void ib_set_iflag();
    void ib_norm_cal(const int cc, const int i, const int j, const int k);
    void nib(const real & n1, const real & n2, const real & n3
           , const real & n4, const real & n5, const real & n6
           , real r[]);
    void wall_norm(const Scalar & g);
    void nwall(const Scalar & g
             , const real & r1, const real & r2, const real & r3
             , const int & i1, const int & i2, const int &i3
             , real r[] );
    void CIPCSLx1(const Scalar & g1, const Scalar & g2);
    void CIPCSLx2(const Scalar & g1, const Vector & v1);
    void CIPCSLx3(const Scalar & g1, const Vector & v1);
    void CIPCSLx4(const Vector & v1, const Scalar & g1);
    void CIPCSLy1(const Scalar & g1, const Scalar & g2);
    void CIPCSLy2(const Scalar & g1, const Vector & v1);
    void CIPCSLy3(const Scalar & g1, const Vector & v1);
    void CIPCSLy4(const Vector & v1, const Scalar & g1);
    void CIPCSLz1(const Scalar & g1, const Scalar & g2);
    void CIPCSLz2(const Scalar & g1, const Vector & v1);
    void CIPCSLz3(const Scalar & g1, const Vector & v1);
    void CIPCSLz4(const Vector & v1, const Scalar & g1);
    void funclxyz(real & r1, real & r2, real r3, real r4, real r5, real r6
                 ,real r7, real r8);
    void bdphiface(const Vector & v, const Comp & m, const Scalar & g1);
    void bdphiface_check(const Vector & v, const Scalar & g1);
    real totalvol(const Scalar & g);
    void plot_f(const char *);
    void plot_sigx(const char *);
    void plot_sigy(const char *);
    void plot_sigz(const char *);
    void plot_sxyz(const char *, const Comp & m);
    void set_minval(real r) {minclr=r;}
    void set_maxval(real r) {maxclr=r;}
    real beta(const real a1, const real a2, const bool b);

    void cal_fs();
    void fs_bnd_nosubgrid(const Scalar & scp);
    real extrapolate_c(const Scalar & sca, const int i, const int j, const int k,
                       const int ofx, const int ofy, const int ofz,
                       const real rat);
    void cal_adens();
    void update_at_walls(Scalar & sca);

    void interfacial_flagging(Scalar & scp);
    bool Interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);
    bool Interface(const int i, const int j, const int k);

    Scalar clr, sclr;            /* color function, smeared color function */

    /* ancillary */
    Heaviside * heavi;
    Vector fs;
    Scalar adens;

    Scalar nx,ny,nz;       /* normal to interface */
    Scalar dist;           /* distance function */
    Scalar kappa;          /* curvature function */
    Scalar alp;            /* coefficient for local sharpening */
    Scalar fn, atmp, stmp; /* temporary */
    ScalarInt iflag, wflag, intflag;
    Scheme scheme;
    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real tol_wall; /* wall tolerance for erroneous interfaces */
    real pi,phimin,phimax,phisurf;
    real dxmin,ww;
    real sum_outlet,sum_outletm,clrsum1,clrsum2;
    int nredist, itsharpen, nlayer, ialpcal, itsmear;
    real *** vel;
    real *** delrho;
    real epsnorm,epss,eps_clr,eps_st;
    Vector sxyz;
    real cangle;
    real minclr, maxclr;
    bool localSharpen, use_dist_for_kappa;
};
#endif
