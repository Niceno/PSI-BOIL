#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#include "update_step.cpp"
#define OPT_RACK
#define PHASECHANGE
//#define RESTART_COPY_TPR
//#define COPY_TPR
using namespace std;

const int  gLevel = 8; //32

/* Domain */
const int  NX1 = 16*gLevel;
const int  NX2 =  2*gLevel;
const int  NY  = 8*gLevel;
const int  NZ  = 8*gLevel;
const int  NmZ =  2*gLevel;
const int  NZA = NZ + NmZ;
const real LY  =  0.005;
const real LX1 =  0.01;
const real LX3 =  0.012;
const real LZ  =  0.00589;
const real LmZ = -0.001;

/* parameter for boundary conditions */
const real prs = 10.5*1e+5;  // system pressure (Pa)
const real tsat0 = 0.0;
const real dtsub = 10.0;
const real qflux = 1.0e+6;  // heater power (W/m2)
const real thickITO = 0.7e-6; // thickness of ITO heater (m)
const real qsrc = qflux/thickITO;  // heater power (W/m3)

/* constants */
const real gravity = 9.8;

double frac(double x) {
  return x - std::floor(x);
}


/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  /*----------+
  |  grid(s)  |
  +----------*/
  const real dx = LX1/real(NX1);
  Grid1D gx1( Range<real>(0.0, LX1), NX1, Periodic::no());
  Grid1D gx2( Range<real>(LX1, LX3), Range<real>(1.2*dx,2*dx), NX2, Periodic::no());
  Grid1D gx( gx1, gx2, Periodic::no());
  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), NY, Periodic::yes());

  Grid1D gz1( Range<real>(LmZ,0.0), Range<real>(-LmZ/NmZ*2.0,thickITO), NmZ, Periodic::no());
  //Grid1D gz2( Range<real>(0.0,LZ), NZ, Periodic::no());
  Grid1D gz2( Range<real>(0.0,0.5*LZ), Range<real>(dx/4.0,dx),NZ*3/4, Periodic::no());
  Grid1D gz3( Range<real>(0.5*LZ,LZ),  Range<real>(1.2*dx,4.0*dx), NZ/4, Periodic::no());
  Grid1D gztmp( gz1, gz2, Periodic::no());
  Grid1D gz( gztmp, gz3, Periodic::no(), BndGrid::wall(), BndGrid::symmetry());

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Domain dom(gx, gy, gz, & floor);
  int nig=gx.ncell();
  int njg=gy.ncell();
  int nkg=gz.ncell();
  //boil::plot->plot(dom,"dom");
  dom.save("grid32.grd");

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(dom), xyz(dom);             // velocity
  Scalar press(dom), p  (dom), f  (dom); // pressure
  Scalar tpr(dom), q  (dom);             // temperature
  Scalar c  (dom), g  (dom), step(dom), kappa(dom); // color function
  Scalar mdot  (dom); // phase change mdot 2
  Scalar wd(dom), mu_t(dom), yplus(dom);// eddy viscosity
  Scalar sdummy(dom);  // temporary or nousage
  Scalar id_bubble(dom);  // for Floodfill
  real *** copy_planeVec;

  /*--------------+
  | fluent-2D.txt |
  +--------------*/
  int nk_fluent;
  std::vector<real> z_fluent;
  std::vector<real> TKE_fluent;
  std::vector<real> len_fluent;
  std::fstream input_TKE;
  input_TKE.open("fluent-2D.txt", std::ios::in);
  if( !input_TKE.fail() ) {
    input_TKE >> nk_fluent;
    boil::oout<<"nk_fluent= "<<nk_fluent<<"\n";
    std::string line;
    input_TKE >> line >> line >> line >> line >> line;
    for(int K=0; K<nk_fluent; K++) {
      real val1,val2,val3,val4,val5;
      input_TKE >> val1 >> val2 >> val3 >> val4 >> val5;
      //boil::oout<< val1<<" "<<val2<<" "<<val5<<"\n";
      z_fluent.push_back(val1);
      TKE_fluent.push_back(val2);
      len_fluent.push_back(val5);
    }
  } else {
    boil::oout<<"cannot find fluent-2D.txt"<<"\n";
    exit(0);
  }
  //for(int K=0; K<nk_fluent; K++) {
  //  boil::oout<<"K,z_fluent,TKE,len"<<K<<" "<<z_fluent[K]<<" "<<TKE_fluent[K]<<" "<<len_fluent[K]<<"\n";
  //}
  // calculate TKE and turbulent length at cell center
  real TKE_local[p.ek()+1], len_local[p.ek()+1];
  for_vk(p,k) {
    //boil::oout<<"k= "<<k<<"\n";
    real ztmp = p.zc(k);
    if (ztmp<=0.0) {
      TKE_local[k]=0.0;
      len_local[k]=0.0;
    } else {
      for(int K=0; K<nk_fluent; K++) {
        if(ztmp<z_fluent[K]) {
          real s1 = ztmp - z_fluent[K-1];
          real s2 = z_fluent[K] - ztmp;
          TKE_local[k] = s2/(s1+s2)*TKE_fluent[K-1]+s1/(s1+s2)*TKE_fluent[K];
          len_local[k] = s2/(s1+s2)*len_fluent[K-1]+s1/(s1+s2)*len_fluent[K];
          //continue;
          break;
        }
      }
    }
  }
  for_vk(p,k) {
    boil::oout<<"k,p.zc,TKE_local,len_local: "<<k<<" "<<p.zc(k)<<" "<<TKE_local[k]<<" "<<len_local[k]<<"\n";
  }
  // calculate TKE and turbulent length at w-cell center
  real TKE_local_w[p.ek()+1], len_local_w[p.ek()+1];
  for_vmk(uvw,Comp::w(),k) {
    //boil::oout<<"k= "<<k<<"\n";
    real ztmp = uvw.zc(Comp::w(),k);
    if (ztmp<=0.0) {
      TKE_local_w[k]=0.0;
      len_local_w[k]=0.0;
    } else {
      for(int K=0; K<nk_fluent; K++) {
        if(ztmp<z_fluent[K]) {
          real s1 = ztmp - z_fluent[K-1];
          real s2 = z_fluent[K] - ztmp;
          TKE_local_w[k] = s2/(s1+s2)*TKE_fluent[K-1]+s1/(s1+s2)*TKE_fluent[K];
          len_local_w[k] = s2/(s1+s2)*len_fluent[K-1]+s1/(s1+s2)*len_fluent[K];
          //continue;
          break;
        }
      }
    }
  }
  for_vk(p,k) {
    boil::oout<<"k,zcw,TKE_local_w,len_local_w: "<<k<<" "<<uvw.zc(Comp::w(),k)<<" "<<TKE_local_w[k]<<" "<<len_local_w[k]<<"\n";
  }
}
