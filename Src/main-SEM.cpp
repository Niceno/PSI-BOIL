#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#include "update_step.cpp"
#include "SEM.cpp"
//#define OPT_RACK
//#define PHASECHANGE
//#define RESTART_COPY_TPR
#define COPY_TPR
using namespace std;

#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

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

real frac(real x) {
  return x - std::floor(x);
}


/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc==1){
    boil::oout<<"One command line argument is required!"<<"\n";
    boil::oout<<"./Boil wmin (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

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
  dom.save("grid8.grd");

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
  real *** copy_planeVec;
  //alloc3d( &copy_planeVec, 3, uvw.nj()+1, uvw.nk()+1);
  //real copy_planeVec[3][uvw.nj()+1][uvw.nk()+1];

#ifdef OPT_RACK
  /*--------------------+
  |  monitoring points  |
  +--------------------*/
  // find I
  real Xr1=0.001, Xr2=0.003, Xr3=0.005, Xr4=0.007, Xr5=0.009, Xr0=dx*0.5;
  int Ir1,Ir2,Ir3,Ir4,Ir5,Ir0;
  for(int I=boil::BW; I<=NX1+boil::BW-1; I++) {
    real xtmp = dom.xc_global(I);
    //boil::oout<<"I "<<I<<" xtmp "<<xtmp<<" "<<dx<<"\n";
    if(Xr1<=xtmp && xtmp<Xr1+dx) Ir1=I-boil::BW;
    if(Xr2<=xtmp && xtmp<Xr2+dx) Ir2=I-boil::BW;
    if(Xr3<=xtmp && xtmp<Xr3+dx) Ir3=I-boil::BW;
    if(Xr4<=xtmp && xtmp<Xr4+dx) Ir4=I-boil::BW;
    if(Xr5<=xtmp && xtmp<Xr5+dx) Ir5=I-boil::BW;
    //if(Xr0<=xtmp && xtmp<Xr0+dx) Ir0=I-boil::BW;
  }
  Ir0=1;
  boil::oout<<"main:Ir1,Ir2,Ir3,Ir4,Ir5,Ir0="<<Ir1<<" "<<Ir2<<" "<<Ir3<<" "<<Ir4<<" "<<Ir5<<" "<<Ir0<<"\n";
  Location m0("m0", dom, Ir1+1, NY/2, NmZ+4);  // 4 cell above wall
  // x~0.001
  Rack r1("r1", dom, Ir1+1, Range<int>(1,NY), Range<int>(1,NZA));
  r1.save_grid(c,"r1-scalar");
  for_vijk(p,i,j,k){p[i][j][k]=p.dV(i,j,k);}
  r1.save(p,"r1-volume",0);
  r1.save_grid(uvw,Comp::u(),"r1-vector-U");
  r1.save_grid(uvw,Comp::v(),"r1-vector-V");
  r1.save_grid(uvw,Comp::w(),"r1-vector-W");

  // x~0.003
  Rack r2("r2", dom, Ir2, Range<int>(1,NY), Range<int>(1,NZA));
  r2.save_grid(c,"r2-scalar");
  r2.save(p,"r2-volume",0);
  r2.save_grid(uvw,Comp::u(),"r2-vector-U");
  r2.save_grid(uvw,Comp::v(),"r2-vector-V");
  r2.save_grid(uvw,Comp::w(),"r2-vector-W");

  // x~0.005
  Rack r3("r3", dom, Ir3, Range<int>(1,NY), Range<int>(1,NZA));
  r3.save_grid(c,"r3-scalar");
  r3.save(p,"r3-volume",0);
  r3.save_grid(uvw,Comp::u(),"r3-vector-U");
  r3.save_grid(uvw,Comp::v(),"r3-vector-V");
  r3.save_grid(uvw,Comp::w(),"r3-vector-W");

  // x~0.007
  Rack r4("r4", dom, Ir4, Range<int>(1,NY), Range<int>(1,NZA));
  r4.save_grid(c,"r4-scalar");
  r4.save(p,"r4-volume",0);
  r4.save_grid(uvw,Comp::u(),"r4-vector-U");
  r4.save_grid(uvw,Comp::v(),"r4-vector-V");
  r4.save_grid(uvw,Comp::w(),"r4-vector-W");

  // x~0.009
  Rack r5("r5", dom, Ir5, Range<int>(1,NY), Range<int>(1,NZA));
  r5.save_grid(c,"r5-scalar");
  r5.save(p,"r5-volume",0);
  p=0.0;
  r5.save_grid(uvw,Comp::u(),"r5-vector-U");
  r5.save_grid(uvw,Comp::v(),"r5-vector-V");
  r5.save_grid(uvw,Comp::w(),"r5-vector-W");

  // x~0.00
  Rack r0("r0", dom, Ir0, Range<int>(1,NY), Range<int>(1,NZA));
  r0.save_grid(c,"r0-scalar");
  r0.save(p,"r0-volume",0);
  p=0.0;
  r0.save_grid(uvw,Comp::u(),"r0-vector-U");
  r0.save_grid(uvw,Comp::v(),"r0-vector-V");
  r0.save_grid(uvw,Comp::w(),"r0-vector-W");
#endif

  /*----------------+
  | fluent-2DR1.txt |
  +----------------*/
  int nk_fluent;
  std::vector<real> z_fluent;
  std::vector<real> U_fluent;
  std::vector<real> tke_fluent;
  std::vector<real> omg_fluent;
  std::vector<real> len_fluent;
  std::vector<real> nut_fluent;
  std::fstream input_tke;
  input_tke.open("fluent-2DR3.txt", std::ios::in);
  if( !input_tke.fail() ) {
    input_tke >> nk_fluent;
    boil::oout<<"nk_fluent= "<<nk_fluent<<"\n";
    std::string line;
    input_tke >> line >> line >> line >> line >> line >> line >> line;
    for(int K=0; K<nk_fluent; K++) {
      real val1,val2,val3,val4,val5,val6,val7;
      input_tke >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
      z_fluent.push_back(val1);
      tke_fluent.push_back(val2);
      omg_fluent.push_back(val3);
      len_fluent.push_back(val5);
      U_fluent.push_back(val6);
      nut_fluent.push_back(val7);
    }
    //exit(0);
  } else {
    boil::oout<<"cannot find fluent-2DR3.txt"<<"\n";
    exit(0);
  }
  //for(int K=0; K<nk_fluent; K++) {
  //  boil::oout<<K<<" "<<z_fluent[K]<<" "<<tke_fluent[K]<<"\n";
  //}

#if 0
  // check how to access global
  boil::oout<<"nkg= "<<nkg<<" dom.gk "<<dom.gk()<<" dom.gik "<<dom.gik()<<"\n";
  // loop for cell
  for (int KC=0; KC<dom.gik(); KC++){
    boil::oout<<"cell: KC= "<<KC<<" dom.zc_global "<<dom.zc_global(KC+boil::BW)<<"\n";
  }
  // loop for node
  for (int KN=0; KN<=dom.gik(); KN++){
    boil::oout<<"node: KN= "<<KN<<" dom.zn_global "<<dom.zn_global(KN+boil::BW)<<"\n";
  }
#endif
  // store U, tke, turbulent length
  real U_global[dom.gik()+1], nut_global[dom.gik()], tke_global[dom.gik()],
       omg_global[dom.gik()], len_global[dom.gik()];
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if (ztmp<=0.0) {
      U_global[KC]=0.0;
      nut_global[KC]=0.0;
      tke_global[KC]=0.0;
      omg_global[KC]=0.0;
      len_global[KC]=0.0;
    } else {
      for(int KF=0; KF<nk_fluent; KF++) {
        if(ztmp<z_fluent[KF]) {
          real s1 = ztmp - z_fluent[KF-1];
          real s2 = z_fluent[KF] - ztmp;
          U_global[KC]   = s2/(s1+s2)*U_fluent[KF-1]  +s1/(s1+s2)*U_fluent[KF];
          nut_global[KC] = s2/(s1+s2)*nut_fluent[KF-1]+s1/(s1+s2)*nut_fluent[KF];
          tke_global[KC] = s2/(s1+s2)*tke_fluent[KF-1]+s1/(s1+s2)*tke_fluent[KF];
          omg_global[KC] = s2/(s1+s2)*omg_fluent[KF-1]+s1/(s1+s2)*omg_fluent[KF];
          len_global[KC] = s2/(s1+s2)*len_fluent[KF-1]+s1/(s1+s2)*len_fluent[KF];
          break;
        }
      }
    }
    //boil::oout<<"KC,Z,U,mut,tke,omg,len "<<KC<<" "<<dom.zc_global(KC+boil::BW)<<" "
    //          <<U_global[KC]<<" "<<nut_global[KC]<<" "
    //          <<tke_global[KC]<<" " <<omg_global[KC]<<" "<<len_global[KC]<<"\n";
  }
  // buffer cell
  U_global[dom.gik()]=U_global[dom.gik()-1];

  // calculate du/dz
  real dudz_global[dom.gik()];
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    real ztmp_minus = dom.zc_global(KC-1+boil::BW);
    if (ztmp<=0.0) {
      dudz_global[KC]=0.0;
    } else if (ztmp>0 && ztmp_minus<0) {
      dudz_global[KC]=(U_global[KC]-0.0)/(ztmp-0.0);
      //boil::oout<<"dudz at wall "<<U_global[KC]<<" "<<ztmp<<" "<<dudz_global[KC]<<" "<<KC<<"\n";
    } else {
      dudz_global[KC]=(U_global[KC+1]-U_global[KC-1])
                     /(dom.zc_global(KC+1+boil::BW)-dom.zc_global(KC-1+boil::BW));
      //if (KC==dom.gik()-1) {
      //  std::cout<<"ZC "<<dom.zc_global(KC+1+boil::BW)<<" "<<dom.zc_global(KC-1+boil::BW)<<"\n";
      //  exit(0);
      //}
    }
    //boil::oout<<"KC,ZC,U,dudz "<<KC<<" "<<dom.zc_global(KC+boil::BW)<<" "
    //          <<U_global[KC]<<" "<<dudz_global[KC]<<"\n";
  }

  // calculate Reynolds stress Rij(KC)
  real L11_global[dom.gik()], L22_global[dom.gik()], L31_global[dom.gik()], L33_global[dom.gik()];
  real sig_global[dom.gik()];
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if (ztmp<=0.0) {
      L11_global[KC]=0.0;
      L22_global[KC]=0.0;
      L31_global[KC]=0.0;
      L33_global[KC]=0.0;
    } else {
      real Rxx,Ryy,Rzz,Rxz,Rxz_max,dz;
      Rxx=2.0/3.0*tke_global[KC];
      Ryy=Rxx;
      Rzz=Rxx;
      Rxz=-nut_global[KC]*dudz_global[KC];
      Rxz_max = 0.98*sqrt(Rxx*Rzz);
      if (fabs(Rxz)>Rxz_max) Rxz=copysign(Rxz_max,Rxz);
      L11_global[KC]=sqrt(Rxx);
      L22_global[KC]=sqrt(Ryy);
      L31_global[KC]=Rxz/sqrt(Rxx);
      L33_global[KC]=sqrt(Rzz-Rxz*Rxz/Rxx);
      dz = dom.zn_global(KC+1+boil::BW)-dom.zn_global(KC+boil::BW);
      sig_global[KC]=max(dz,len_global[KC]);
      //boil::oout<<"KC,dz,len,sig "<<KC<<" "<<dz<<" "<<len_global[KC]<<" "<<sig_global[KC]<<"\n";
    }
  }

  // representative sig_global
  real sum_sig=0.0;
  int n_sig=0;
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if(ztmp>0.0){
      sum_sig+=sig_global[KC];
      n_sig++;
    }
  }
  const real sig_rep = sum_sig/real(n_sig);

  // Number of eddies N_eddy
  int N_eddy = int(4.0 * (LY*LZ)/(sig_rep*sig_rep));
  boil::oout<<"N_eddy= "<<N_eddy<<" sig_rep "<<sig_rep<<"\n";

  // ----- SEM state (global-ish in main) -----
  std::mt19937_64 rng(1234567ULL);         // or seed from time
  std::vector<Eddy> E;                     // eddy list
  BoundsYZ B{ -LY/2.0, +LY/2.0, 0.0, LZ }; // your inlet plane extents (y,z)

  // initialize eddies uniformly in source slab [-sigma_rep, 0] × [ymin,ymax] × [zmin,zmax]
  E.resize(N_eddy);
  for (int n=0; n<N_eddy; ++n){
    E[n].x  = urand(rng, -sig_rep, 0.0);
    E[n].y0 = urand(rng, B.y_min, B.y_max);
    E[n].z0 = urand(rng, B.z_min, B.z_max);
    for (int i=0;i<3;++i) E[n].sgn[i] = rand_sign(rng);
  }
  for (int i=0; i<N_eddy; i+=100) {
    boil::oout<<"E[i]= "<<i<<" "<<E[i].x<<" "<<E[i].y0<<" "<<E[i].z0<<"\n"; 
  }

  const real Uc = 1.127234319;
  boil::oout<<"dt recommended based on Uc "<<0.2*sig_global[NmZ]/Uc<<" "<<NmZ<<"\n";
  // Precompute a fixed “alpha” for EWMA (dt will be adapted later)
  //const real Tref = sig_rep / Uc;              // eddy lifetime
  //const real tau  = 8.0 * Tref;                // averaging horizon ~ 510 T
  //real alpha_EWMA = (1e-5) / tau;              // will update with actual dt later
  //boil::oout<<"alpha_EWMA= "<<alpha_EWMA<<"\n";
  //exit(0);

  // Allocate EWMA accumulators per inlet node (component × j × k)
  std::vector<EWMA> ew_u((uvw.nj()+1)*(uvw.nk()+1));
  std::vector<EWMA> ew_v((uvw.nj()+1)*(uvw.nk()+1));
  std::vector<EWMA> ew_w((uvw.nj()+1)*(uvw.nk()+1));
  auto idx = [&](int j,int k){ return j*(uvw.nk()+1)+k; };

  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    //uvw.bc(m).add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell()),
    //          Range<int>(1,NmZ), BndType::inlet(), 0.0, 0.0, 0.0 ) );
    //uvw.bc(m).add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell()), Range<int>(NmZ+1,gz.ncell()), 
    //          BndType::inlet(), "1.3201*(1.0-((0.00589-z)/0.00589)^1.2)^0.2", 0.0, 0.0 ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::insert() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  }
  press.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  f    = p.shape();
  q    = p.shape();
  mdot = p.shape();
  kappa = p.shape();
  g     = p.shape();
  yplus = p.shape();
  mu_t  = p.shape();
  wd    = p.shape();
  sdummy  = p.shape();

  tpr.bc().add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell())
            , Range<int>(1,NmZ), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell())
            , Range<int>(NmZ+1,gz.ncell()), BndType::dirichlet(), tsat0-dtsub) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );

  c.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 1.0 ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  step = c.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(dom), liquid(dom), sapphire(dom); // zircaloy(dom)
  //zircaloy.lambda  (13.94);
  //zircaloy.rho    (6537.);
  //zircaloy.cp     (295.9*6537.);
  sapphire.rho    (3980.0);
  sapphire.cp     (750.0*3980.0);
  sapphire.lambda (35.0);

  const real tsat = IF97::Tsat97(prs);
  vapor.mu(IF97::viscvap_p(prs));
  vapor.rho(IF97::rhovap_p(prs));
  vapor.cp(IF97::cpvap_p(prs)*IF97::rhovap_p(prs));
  vapor.lambda(IF97::tcondvap_p(prs));

  liquid.mu(IF97::viscliq_p(prs));
  liquid.rho(IF97::rholiq_p(prs));
  liquid.cp(IF97::cpliq_p(prs)*IF97::rholiq_p(prs));
  liquid.lambda(IF97::tcondliq_p(prs));
  Matter mixed(liquid, vapor, &step);
  mixed.sigma(IF97::sigma97(tsat));
  const real latent=IF97::hvap_p(prs)-IF97::hliq_p(prs);

  /* estimate betal */
  const real liquid_drhodt = (IF97::rholiq_p(1.01*prs)-IF97::rholiq_p(0.99*prs))
                             /(IF97::Tsat97(1.01*prs)-IF97::Tsat97(0.99*prs));
  const real vapor_drhodt=0.0; //[kg/m3K]

  boil::oout<<"properties at pressure "<<prs<<boil::endl;
  boil::oout<<"Tsat (K)= "<<tsat<<" vapor= "<<IF97::viscvap_p(prs)<<" (Pa.s) "
           <<IF97::rhovap_p(prs)<<" (kg/m3) "<<IF97::cpvap_p(prs)<<" (J/kg.K) "
	   <<IF97::tcondvap_p(prs)<<" (W/m.K)\n";
  boil::oout<<"Tsat (C)= "<<tsat-273.15<<" liquid= "<<IF97::viscliq_p(prs)<<" (Pa.s) "
           <<IF97::rholiq_p(prs)<<" (kg/m3) "<<IF97::cpliq_p(prs)<<" (J/kg.K) "
           <<IF97::tcondliq_p(prs)<<" (W/m.K)\n";
  boil::oout<<"sigma= "<<IF97::sigma97(tsat)<<" latent= "<<latent
           <<" drho/dT= "<<liquid_drhodt<<"\n";

  const real rhol  = liquid.rho()->value();
  const real rhov  = vapor.rho()->value();
  const real mul   = liquid.mu()->value();
  const real muv   = vapor.mu()->value();
  const real sigma = mixed.sigma()->value();
  real const mu_t_max =100.0*std::max(mul,muv);  // previously 10.0

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 50000000;
  const real tint = 1.0e-3;  // time interval for dat-output
  const real tint2 = 2.0e-5;  //2.0e-5; // time interval for monitoring rack
  const int  nint = 10000;   // time interval for bck-output
  const real dxmin = dom.dxyz_min(); // minimum size of cell in liquid
#ifdef PHASECHANGE
  const real dt  = 10.0*pow(vapor.rho()->value()*pow(dxmin,3.0)
	             / (2.0*3.1415*sigma),0.5);
#else
  const real dt = 1e-3;
#endif
  const real dt0 = 1e-5;
  const real cfl_limit=0.25;
  Times time(ndt, dt0);
  time.print_time(false);
  time.set_coef_dec(0.5);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(dom, Prec::ic2());
  Krylov * solver2 = new BiCGS(dom, Prec::di());

  Pressure pr( p,   f,   uvw, time, solver, &mixed );
  Momentum ns( uvw, xyz,        time, solver, &mixed );
  ns.convection_set(TimeScheme::adams_bashforth()); // Grid32R13-woPC
  ns.convection_set(ConvScheme::central()); // Grid32R12-woPC
  ns.diffusion_set(TimeScheme::crank_nicolson()); // Grid32R13-woPC
  EnthalpyFD en(tpr, q, c, uvw, time, solver2, & mixed, tsat0, & sapphire);
  en.convection_set(TimeScheme::forward_euler());
  en.diffusion_set(TimeScheme::backward_euler());
  CIPCSL2 conc (c,  g, kappa, uvw, time, solver);
  conc.set_itsharpen(10);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_cangle(0.0);

  const real rseed = 1.2 * dx;
  const real dmicro_min = 1.0e-10; // nonuse
#ifdef PHASECHANGE
  Nucleation nucl( &c, &tpr, &q, &time, sdummy, &mixed,
    rseed, dmicro_min, latent, conc.get_cangle(), tsat0, &sapphire);
  nucl.set_seed_period(1e-5);    // continue to plant color function
  nucl.set_pre_heat_sink(true);  // give heat sink before plant of color function
  nucl.set_micro_exists(false);  // micro-layer model is not used
  nucl.set_rseed_plus(2.0);      // clr detect range for replant

  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw,
                 time, &mixed, latent, tsat0, &sapphire, &nucl);
#endif

  AC sol( &pr );
  sol.stop_if_diverging(true);
  sol.min_cycles(3);
  sol.max_cycles(10);

  /* wall distance */
  Distance di(wd, sdummy, uvw, time, solver);

  /* eddy viscosity */
  Model tm;

  /*-------------------+
  |  check if restart  |
  +-------------------*/

  std::fstream input;
  int irun = 0;
  if(boil::cart.iam()==0){
    input.open("run.txt", std::ios::in);
    if( !input.fail() ) {
      input >> irun;
      boil::oout<<"read irun.  irun= "<<irun<<"\n";
    }
    input.close();
  }
  boil::cart.sum_int(&irun);
  if (irun==1){
    boil::oout<<"exit job due to irun=1"<<"\n";
    exit(0);
  }

  if(boil::cart.iam()==0){
    std::fstream output;
    output.open("run.txt", std::ios::out);
    output << 1 << boil::endl;
    output.close();
  }

  int ts=0;
  bool restart = false;

  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    restart=true;
  }

  if( restart ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
  }
 
  if( restart ) {
    uvw   .load("uvw", ts);
    press .load("press", ts);
    tpr   .load("tpr", ts);
#if 1
    /* normal start */
    conc  .load("conc", ts);
#ifdef PHASECHANGE
    nucl  .load("nucl", ts);
#endif
    wd    .load("wd",   ts);
#else
    /* after changeGrid */
    /* load conc-phi as a scalar */ 
    c.load("conc-phi", ts);
    c.exchange_all();
    /* calculate conc-sigma etc. */ 
    conc.init();
    /* zsite can be changed */
#ifdef PHASECHANGE
    real zsite=rseed*cos(conc.get_cangle()/180.0*boil::pi);
    boil::oout<<"main:change_site_z= "<<zsite<<" "<<rseed<<"\n";
    nucl  .load("nucl", ts, &zsite);
    /* distance function from wall */
#endif
    wd    .load("wd",   ts);  // load as an initial value
    sdummy = 0.0;             // reset source term
    di.compute();
#endif
    /* calculate step in case of restart */
    update_step(c, step, sdummy);
#ifdef RESTART_COPY_TPR
    boil::oout<<"main:tpr.sk() "<<tpr.sk()<<" NmZ "<<NmZ<<"\n";
    for_vijk(tpr,i,j,k){
      if(tpr.xc(i)>0.0006 && tpr.zc(k)<0.0){
        tpr[i][j][k]=tpr[i][j][tpr.sk()+NmZ-1];
      }
    }
    tpr.bnd_update();
    tpr.exchange();
    boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,yplus
                    ,"uvw-c-tpr-press-mdot-mut-yplus",999999);
#endif
  }

  else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    /*--------------------+
    |  initial condition  |
    +--------------------*/
    /* compute distance function */
    sdummy = 0.0; // reset source term for di
    di.compute();
    //boil::plot->plot(wd,"wd",0);

    /* velocity */
    Comp m = Comp::u();
    for_vmijk(uvw,m,i,j,k){
      real zz = uvw.zc(m,k);
        if(zz<0.0) {
          uvw[m][i][j][k] = 0.0;
        } else {
          //uvw[m][i][j][k] = 1.3201*(pow(1.0-pow(((0.00589-zz)/0.00589),1.2),0.2));
          int K = dom.global_K(k)-boil::BW;
          uvw[m][i][j][k] = U_global[K];
        }
    }

    /* color function */
    for_vijk(c,i,j,k)
      c[i][j][k] = 1.0;
    c.exchange_all();
    update_step(c, step, sdummy);
    conc.init();

    /* temperature */
    tpr = tsat0-dtsub;
    for_vijk(tpr,i,j,k){
      if(tpr.zc(k)<0.0){
        tpr[i][j][k] = tsat0;
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();

    /* set nucleation sites */
#ifdef PHASECHANGE
#if 1
    ifstream fin("site-G1000.txt");
    if( fin.fail() ) {
      boil::oout<<"Cannot open nucleation site file\n";
      exit(0); 
    } else {
      boil::oout<<"Open nucleation site file\n";
      string ss; 
      getline(fin, ss); // skip first line
      real nsx,nsy,Tact;
      real zsite=rseed*cos(conc.get_cangle()/180.0*boil::pi);
      boil::oout<<"main:zsite= "<<zsite<<" rseed= "<<rseed<<"\n";
      int ins=0;
      std::string line;
      while (std::getline(fin, line)){
        std::istringstream iss(line);
        real nsx, nsy, Tact;
        iss >> nsx >> nsy >> Tact;
        nucl.add(Site(nsx, nsy, zsite, Tact, -0.002));
      	// -0.002 is nonuse: bottom of previous bubble
        if(ins%1000==0){
          boil::oout<<"main:nucl.add= "<<ins<<" "<<nsx<<" "<<nsy<<" "<<Tact<<"\n";
        }
        ins++;
      }
    }
    fin.close();
#else
    #include "nsd.cpp"
#endif
#endif
    /* plot initial condition */
    boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,yplus
                    ,"uvw-c-tpr-press-mdot-mut-yplus",0);
  }
  input.close();

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;
  int iint2 = int(time.current_time()/tint2) + 1;

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# DT:        " << time.dt() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "########################" << boil::endl;

    /*-----------------+
    |  inlet velocity  |
    +-----------------*/
    boil::oout<<"main:inlet_velocity:Hello "<<boil::cart.iam()<<"\n";
    alloc3d( &copy_planeVec, 3, uvw.nj()+1, uvw.nk()+1);
    //real tt = time.current_time();
    real alpha_EWMA = std::min(0.2, time.dt() / (8.0 * (sig_rep/Uc))); // keep α ≤ 0.2
    //1) Convect & recycle eddies
    for (auto & e : E){
      e.x += Uc * dt;
      if (e.x > 0.0){
        e.x  = -sig_rep;
        e.y0 = urand(rng, B.y_min, B.y_max);
        e.z0 = urand(rng, B.z_min, B.z_max);
        for (int i=0;i<3;++i) e.sgn[i] = rand_sign(rng);
      }
    }

    for_m(m){
      //boil::oout<<"main:inlet_velocity:m= "<<m<<" "<<boil::cart.iam()<<"\n";
      //boil::oout<<"main:inlet_velocity:m= "<<boil::cart.iam()<<"\n";
      for_vmjk(uvw,m,j,k) {
        const real zz = uvw.zc(m,k);
        //if(time.current_step()==2)
        //std::cout<<"m,j,k "<<boil::cart.iam()<<" "<<m<<" "<<j<<" "<<k<<" "<<zz<<"\n";
        if(zz<=0.0) {  // solid
          copy_planeVec[~m][j][k]=0.0;
        } else {
#if 0
          real yy = uvw.yc(m,j);
          //real u_mean = 1.3201*(pow(1.0-pow(((0.00589-zz)/0.00589),1.2),0.2));
          //real u_mean = U_local[k];
          //std::cout<<" u_mean "<<k<<" "<<u_mean<<" "<<U_local[k]<<"\n";
          if(m==Comp::u()) {
            real u_mean = U_local[k];
            real amp = sqrt(2.0/3.0*tke_local[k]);
            real omg = 2.0*boil::pi*u_mean/len_local[k];
            real alp = 2.0*boil::pi/len_local[k];
            real phase_diff = 2*boil::pi*frac(sin(alp*yy+alp*zz));
            real u_ = amp*sin(omg*tt+phase_diff);
            //copy_planeVec[~m][j][k] = u_mean;
            copy_planeVec[~m][j][k] = u_mean + u_;
          } else if (m==Comp::v()) {
            real u_mean = U_local[k];
            real amp = sqrt(2.0/3.0*tke_local[k]);
            real omg = 2.0*boil::pi*u_mean/len_local[k];
            real alp = 2.0*boil::pi/len_local[k];
            real dlt = alp/20.0;
            real phase_diff = 2*boil::pi*frac(sin((alp+dlt)*yy+(alp+dlt)*zz));
            real v_ = amp*sin(omg*tt+phase_diff);
            //copy_planeVec[~m][j][k] = 0.0;
            copy_planeVec[~m][j][k] = v_;
          } else if (m==Comp::w()) {
            real u_mean = U_local_w[k];
            real amp = sqrt(2.0/3.0*tke_local_w[k]);
            real omg = 2.0*boil::pi*u_mean/len_local_w[k];
            real alp = 2.0*boil::pi/len_local_w[k];
            real dlt = alp/20.0;
            real phase_diff = 2*boil::pi*frac(sin((alp-dlt)*yy+(alp-dlt)*zz));
            real w_ = amp*sin(omg*tt+phase_diff);
            //copy_planeVec[~m][j][k] = 0.0;
            copy_planeVec[~m][j][k] = w_;
          }
#endif
#if 1
          const real yy = uvw.yc(m,j);

          // local σ at z-index KC (global index: k + BW offset ⇒ use KC= k-? if needed)
          const int KC = dom.global_K(k)-boil::BW;
          const real sigma_loc = sig_global[KC];

          // raw sums s_i(y,z)
          real s_raw[3] = {0.0, 0.0, 0.0};
          const real inv_sqrtN = 1.0 / std::sqrt(real(E.size()));

          for (const auto &e : E){
            const real dx = -e.x;                // plane is at x=0
            const real dy = yy - e.y0;
            const real dz = zz - e.z0;
            const real kern = gauss_iso(dx, dy, dz, sigma_loc);

            s_raw[0] += (real)e.sgn[0] * kern;
            s_raw[1] += (real)e.sgn[1] * kern;
            s_raw[2] += (real)e.sgn[2] * kern;
          }
          s_raw[0] *= inv_sqrtN;
          s_raw[1] *= inv_sqrtN;
          s_raw[2] *= inv_sqrtN;

          // EWMA normalization per-node/component → η
          const int id = idx(j,k);
          real eta_u = ew_u[id].update(s_raw[0], alpha_EWMA);
          real eta_v = ew_v[id].update(s_raw[1], alpha_EWMA);
          real eta_w = ew_w[id].update(s_raw[2], alpha_EWMA);

          // Map to target stresses via Cholesky (axes: x,u  y,v  z,w)
          const real L11 = L11_global[KC];
          const real L22 = L22_global[KC];
          const real L31 = L31_global[KC];
          const real L33 = L33_global[KC];

          const real uxp = L11*eta_u;
          const real uyp = L22*eta_v;
          const real uzp = L31*eta_u + L33*eta_w;

          // Final inlet: mean + fluctuation
          if (m==Comp::u())      copy_planeVec[~m][j][k] = U_global[KC] + uxp;
          else if (m==Comp::v()) copy_planeVec[~m][j][k] =                uyp;
          else                   copy_planeVec[~m][j][k] =                uzp;
#endif
        }
      }
    }
    // copy_planeVec is computed in all processes but copied to uvw only at inlet
    uvw.bnd_insert(Dir::imin(),copy_planeVec);
    //boil::plot->plot(uvw,tpr,"uvw-tpr",1);

    /*---------------------+
    |  reset source terms  |
    +---------------------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /* heat source */
    for_vijk(q,i,j,k)
      q[i][j][k] = 0.0;

    for_vk(q,k){
      if(-thickITO<q.zc(k) && q.zc(k)<0.0){
        for_vij(q,i,j){
          q[i][j][k] = qsrc*q.dV(i,j,k);  // [W/m3]*[m3]=[W]
        }
      }
    }

    /*------------------+
    |  turbulence model |
    +------------------*/
    /* Smagorinsky model */
    //tm.smagorinsky( & ns, & mu_t, 0.17);  // without damping
    const real dist_max = 1.0e-3; // at z = 1mm, yplus ~ 180. damping ~ 0 if yplus > 40.
    tm.smagorinsky( & ns, & mu_t, 0.17, &dist_max, &wd, &yplus);  // with damping
    /* WALE model */
    //tm.wale( & ns, & mu_t, 0.1);  // 0.325^2 ~ 0.544^2
    /* viscous term in wall adjacent cells */
    tm.tau_wall( & ns, wd, & xyz );

    /* limit eddy viscosity */
    for_avijk(mu_t,i,j,k) {
      mu_t[i][j][k] = std::min(mu_t[i][j][k],mu_t_max);
    }

    /*---------------+
    |  phase change  |
    +---------------*/
#ifdef PHASECHANGE
    pc.update( &mu_t );
    //pc.micro(&xyz);
    update_step(c, step, sdummy);  // 0119 need to update step for with IB
    ns.vol_phase_change(&f);
#endif

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* gravity force */
    Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k){
      real phil=step[i][j][k];
      real phiv=1.0-phil;
      real deltmp = tpr[i][j][k]-tsat0;
      real rhomix = (rhol + liquid_drhodt*deltmp)*phil
                  + (rhov + vapor_drhodt*deltmp)*phiv;
      if(dom.ibody().on(m,i,j,k))
        xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
    }
    xyz.exchange();

#ifdef PHASECHANGE
    /* surface tension */
    conc.tension(&xyz, mixed, step);
#endif

#if 1
    /* increase viscosity in outlet region */  // COPY
    const real x0 = LX3 - 0.002;
    const real x1 = LX3 - 0.000;
    for_avi(c,i){
      if(c.xc(i)>x0){
        real coef=std::min((c.xc(i)-x0)/(x1-x0),1.0);
        for_avjk(c,j,k){
          mu_t[i][j][k] = coef * liquid.mu()->value() * 10;
        }
      }
    }
#endif
     
    /* intermediate velocity */
    ns.discretize( &mu_t );

    pr.discretize();
    pr.coarsen();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1.0e-32));

    /* pressure poisson equation */
    p = 0.0;
    sol.vcycle(ResRat(1.0e-16));

    /* update velocity and pressure */
    ns.project(p);
    press += p;

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(dom.ibody().on(i,j,k)){
        if(pmin>press[i][j][k]) pmin=press[i][j][k];
      }
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k){
      if(dom.ibody().on(i,j,k)){
        press[i][j][k] -= pmin;
      } else {
        press[i][j][k] = 0.0;
      }
    }
    press.bnd_update();
    press.exchange_all();

#if 1
    /* limit velocity */
    for_m(m)
      for_avmijk(uvw,m,i,j,k) {
        uvw[m][i][j][k] = max(-25.0,min(25.0,uvw[m][i][j][k]));
      }
#endif

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
#ifdef PHASECHANGE
    conc.advance();
    conc.totalvol();
#endif

    /*---------------------------+
    |  replant seed or cut neck  |
    +---------------------------*/ 
#ifdef PHASECHANGE
    for(int nsd=0; nsd<nucl.size(); nsd++){
      nucl.sites[nsd].set_allow_replant(true);
    }
    nucl.replant();
#endif

    /*------------------------------+
    |  outlet region: delete liquid |
    +------------------------------*/
#ifdef PHASECHANGE
    const real x0b = LX3 - 0.001;
    const real x1b = LX3;
    for_avi(c,i){
      if(c.xc(i)>x0b){
        real coef=std::min((c.xc(i)-x0b)/(x1b-x0b),1.0);
        for_avjk(c,j,k){
          c[i][j][k]= (1.0-coef)*c[i][j][k] + coef*0.0;
        }
      }
    }
    /* update clr in cipcsl2 after seed, cutneck and outlet-region */
    c.bnd_update();
    c.exchange_all();
    conc.update_node(c);
    update_step(c, step, sdummy);

    /* output min & max of color function */
    conc.color_minmax();
    boil::oout<<"main:color_min,max= "<<time.current_time()<<" "
             <<conc.minval()<<" "<<conc.maxval()<<"\n";
#endif

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    en.discretize( &mu_t );
    en.new_time_step( &mu_t );
    en.solve(ResRat(1e-16),"enth");
    //enthFD.solve_sor(24, 1.2, "enthFD");

    /* outlet region */
    for_vi(tpr,i){
      if(c.xc(i)>LX1){
          for_vk(tpr,k) {
            if (tpr.zc(k)<0) {
              for_vj(tpr,j) {
                tpr[i][j][k]=tpr[i-1][j][k];
              }
            } else {
              for_vj(tpr,j) {
                tpr[i][j][k]=max(tpr[i-1][j][k],tsat0);
              }
            }
         }
      }
    }

#ifdef COPY_TPR
    //boil::oout<<"main:tpr.sk() "<<tpr.sk()<<" NmZ "<<NmZ<<"\n";
    for_vijk(tpr,i,j,k){
      if(tpr.xc(i)>0.0006 && tpr.zc(k)<0.0){
        tpr[i][j][k]=tpr[i][j][tpr.sk()+NmZ-1];
      }
    }
    tpr.bnd_update();
    tpr.exchange();
#endif

    tpr.bnd_update();
    tpr.exchange_all();

    /*-------------+
    |  dt control  |
    +-------------*/
    real cflmax = ns.cfl_max();
    boil::oout<<"main:cflmax= "<<time.current_time()<<" "<<cflmax<<"\n";
    /* interface is always included because of outlet */
    time.control_dt(cflmax, cfl_limit, dt);

    if (time.dt()<1e-8) {
      boil::oout<<"Too small time step: "<<time.dt()<<"\n";
      exit(0);
    }

    /*--------------+
    |  output data  |
    +--------------*/
    /* heat flux and area */
    en.hflux_wall_ib( Range<real>(0.0, LX1), Range<real>(-LY/2.0, LY/2.0),&mu_t);

    /* average temperature */
    real sum_tc=0.0;
    real sum_tcb=0.0;
    real vol_tcb=0.0;
    int n_tc=0;
    for_vk(c,k) {
      if(-thickITO<c.zc(k) && c.zc(k)<0.0){
        for_vij(c,i,j) {
          if(c.xc(i)<LX3-0.002) {
            sum_tcb+=tpr[i][j][k]*tpr.dV(i,j,k);
            vol_tcb+=tpr.dV(i,j,k);
          }
        }
      }
    }
    boil::cart.sum_real(&sum_tcb);
    boil::cart.sum_real(&vol_tcb);
    int ktc = c.akm(-1e-7, boil::pico);  // cell include z=-1e-6
    //int ktc = c.akm(LmZ+1e-6, boil::pico);  // cell include z=LmZ+1e-6
    if(ktc>=c.sk() && ktc<=c.ek()) {
      for_vij(c,i,j) {
        if(c.xc(i)<LX3-0.001) {
          sum_tc+=tpr[i][j][ktc];
          n_tc++;
        }
      }
    }
    boil::cart.sum_real(&sum_tc);
    boil::cart.sum_int(&n_tc);

#if 1
    real tc[7];
    for(int ii=1; ii<=6; ii++) {
      real x0=0.0;
      real x_tc;
      x_tc = x0 + (LX3-0.001)/7.0*real(ii);
      if(time.current_step()==1){
        boil::oout<<"main:x_tc= "<<ii<<" "<<x_tc<<"\n";
      }
#if 0
      if (ii==1) { x_tc = x0 + 0.03; }
      else if (ii==2) { x_tc = x0 + 0.06; }
      else if (ii==3) { x_tc = x0 + 0.09; }
      else if (ii==4) { x_tc = x0 + 0.12; }
      else if (ii==5) { x_tc = x0 + 0.15; }
      else { x_tc = x0 + 0.317; }
#endif
      real sum_tc1=0.0;
      int n_tc1=0;
      int ktc = c.akm(-1e-7, boil::pico);  // cell include z=-1e-7
      //int ktc = c.akm(LmZ+1e-6, boil::pico);  // cell include z=LmZ+1e-6
      //std::cout<<"ktc= "<<ktc<<" "<<tpr.zc(ktc)<<" "<<tpr.zc(ktc+1)<<"\n";
      for_vi(tpr,i) {
        if(tpr.xc(i)>=x_tc && tpr.xc(i)<=x_tc+dx) {
          if(ktc>=c.sk() && ktc<=c.ek()) {
            for_vj(tpr,j) {
              sum_tc1+=tpr[i][j][ktc];
              n_tc1++;
            }
          }
        }
      }
      boil::cart.sum_real(&sum_tc1);
      boil::cart.sum_int(&n_tc1);
      tc[ii]=sum_tc1/real(n_tc1);
    }
    std::cout.setf(std::ios_base::scientific);
    std::cout<< std::setprecision(12);
    boil::oout<<"twall= "<<time.current_time()<<" "<<sum_tcb/vol_tcb<<" "
              <<sum_tc/real(n_tc)<<" "<<tc[1]<<" "<<tc[2]<<" "
              <<tc[3]<<" "<<tc[4]<<" "<<tc[5]<<" "<<tc[6]<<"\n";
    std::cout<< std::setprecision(8);
    std::cout.unsetf(std::ios_base::floatfield);
#endif

    /* Sum(mdot) in wall adjascent cells = contact line evaporation */
#ifdef PHASECHANGE
    real sum_mdot_wall_pos = 0.0;
    real sum_mdot_wall_neg = 0.0;
    for_vk(mdot,k) {
      if( approx(mdot.zn(k),0.0,1e-10)) {
        for_vij(mdot,i,j) {
          if(mdot[i][j][k]>0.0) {
            sum_mdot_wall_pos += mdot[i][j][k]*mdot.dV(i,j,k);
          } else {
            sum_mdot_wall_neg += mdot[i][j][k]*mdot.dV(i,j,k);
          }
        }
      }
    }
    boil::cart.sum_real(&sum_mdot_wall_pos);
    boil::cart.sum_real(&sum_mdot_wall_neg);
    boil::oout<<"sum_mdot_wall= "<<time.current_time()<<" "
              <<sum_mdot_wall_pos<<" "<<sum_mdot_wall_neg<<"\n";
#endif

#ifdef OPT_RACK
    /* monitoring point, every time step */
    boil::oout<<"m0.val= "<<time.current_time()<<" "<<time.dt()<<" c= "
              <<m0.value(c)<<" tpr= "<<m0.value(tpr)<<" U= "<<m0.value(uvw,Comp::u())<<" "
              <<m0.value(uvw,Comp::v())<<" "<<m0.value(uvw,Comp::w())<<"\n";
#endif

    /* tecplot files */
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      tpr.exchange_all();
      boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,yplus,
                      "uvw-c-tpr-press-mdot-mut-yplus",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /* monitoring */
#ifdef OPT_RACK
    boil::oout<<"main:iint2= "<<int(time.current_time()/(tint2))<<" "<<iint2<<"\n";
    if(int(time.current_time()/(tint2)) >= iint2 ) {
      iint2 = int(time.current_time() / tint2);
      tpr.exchange_all();
      boil::oout<<"main:Plotting:Rack\n";
      // r1
      r1.save(uvw,Comp::u(),"r1-U",iint2);
      r1.save(uvw,Comp::v(),"r1-V",iint2);
      r1.save(uvw,Comp::w(),"r1-W",iint2);
      r1.save(c,"r1-c",iint2);
      r1.save(tpr,"r1-tpr",iint2);
      r1.save(mu_t,"r1-mut",iint2);
      // r2
      r2.save(uvw,Comp::u(),"r2-U",iint2);
      r2.save(uvw,Comp::v(),"r2-V",iint2);
      r2.save(uvw,Comp::w(),"r2-W",iint2);
      r2.save(c,"r2-c",iint2);
      r2.save(tpr,"r2-tpr",iint2);
      r2.save(mu_t,"r2-mut",iint2);
      // r3
      r3.save(uvw,Comp::u(),"r3-U",iint2);
      r3.save(uvw,Comp::v(),"r3-V",iint2);
      r3.save(uvw,Comp::w(),"r3-W",iint2);
      r3.save(c,"r3-c",iint2);
      r3.save(tpr,"r3-tpr",iint2);
      r3.save(mu_t,"r3-mut",iint2);
      // r4
      r4.save(uvw,Comp::u(),"r4-U",iint2);
      r4.save(uvw,Comp::v(),"r4-V",iint2);
      r4.save(uvw,Comp::w(),"r4-W",iint2);
      r4.save(c,"r4-c",iint2);
      r4.save(tpr,"r4-tpr",iint2);
      r4.save(mu_t,"r4-mut",iint2);
      // r5
      r5.save(uvw,Comp::u(),"r5-U",iint2);
      r5.save(uvw,Comp::v(),"r5-V",iint2);
      r5.save(uvw,Comp::w(),"r5-W",iint2);
      r5.save(c,"r5-c",iint2);
      r5.save(tpr,"r5-tpr",iint2);
      r5.save(mu_t,"r5-mut",iint2);
      // increase counter
      // r0
      r0.save(uvw,Comp::u(),"r0-U",iint2);
      r0.save(uvw,Comp::v(),"r0-V",iint2);
      r0.save(uvw,Comp::w(),"r0-W",iint2);
      r0.save(c,"r0-c",iint2);
      r0.save(tpr,"r0-tpr",iint2);
      r0.save(mu_t,"r0-mut",iint2);
      // increase counter
      iint2 = int(time.current_time()/tint2) + 1;
    }
#endif

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % nint==0) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
#ifdef PHASECHANGE
      nucl .save("nucl",   time.current_step());
#endif
      wd   .save("wd",   time.current_step());
    }


    if( boil::timer.current_min() > (wmin)
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
#ifdef PHASECHANGE
      nucl .save("nucl",   time.current_step());
#endif
      wd   .save("wd",   time.current_step());

      uvw  .rm("uvw",   ts);
      press.rm("press", ts);
      conc .rm("conc",  ts);
      tpr  .rm("tpr",   ts);
#ifdef PHASECHANGE
      nucl .rm("nucl", ts);
#endif
      wd   .rm("wd", ts);
    }

    
    if((time.current_step()) % (nint)==0 ) {
      if( boil::cart.iam()==0) {
        std::fstream output;
        std::stringstream ss;
        ss <<"time-"<<time.current_step()<<".txt";
        std::string fname = ss.str();
        int len = fname.length();
        char * cfname = new char[len+1];
        memcpy(cfname, fname.c_str(), len+1);
        output << std::setprecision(16);
        output.open(cfname, std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output.close();
      }
    }
 
    if( boil::timer.current_min() > (wmin)
      || time.current_step()==time.total_steps()) {
      std::fstream output;
      output << std::setprecision(16);
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      output.open("run.txt", std::ios::out);
      output << 0 << boil::endl;
      output.close();
      boil::timer.stop();
      boil::timer.report();
      exit(0); 
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}
/*-----------------------------------------------------------------------------+
 '$Id: main-dam.cpp,v 1.12 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
