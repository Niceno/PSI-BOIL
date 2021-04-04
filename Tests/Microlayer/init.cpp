  int ts,iint(0);
  /* load variables */
  std::vector<Scalar*> load_scalars = { &press, &c, &tpr };
  std::vector<std::string> load_scalar_names = { "press", "c", "tpr" };

  std::vector<Vector*> load_vectors = { &uvw };
  std::vector<std::string> load_vector_names = { "uvw" };

  std::vector<Nucleation*> load_nucls = {};
  std::vector<std::string> load_nucl_names = {};
  std::vector<CIPCSL2*> load_cipcsls = {};
  std::vector<std::string> load_cipcsl_names = {};

  std::vector<real*> store_reals = {};
  std::vector<int*>  store_ints  = { &iint };

  const real heater_extent = 1.81e-3;
  const real outer_factor = 0.;

  auto heatfunc = [&](const real x, const real y, const real c, const real v) {
    real h = qsrc*c*v;
    if(case_flag>0)
      return h;

#if 0
    if(fabs(x)<=heater_extent&&fabs(y)<=heater_extent) {
    } else {
      h *= outer_factor;
    }
#else
    h *= std::max(0.,(4.0e-3-x)/4.0e-3);
#endif
  
    return h;
  };

  /* solid */
  csub = 0.0;
  for_vijk(csub,i,j,k) {
    if(csub.zc(k)<0.) {
      if(csub.zn(k+1)<=-LZheat) {
        csub[i][j][k] = 1.0;
      } else if(csub.zn(k)<= -LZheat) {
        csub[i][j][k] = (fabs(csub.zn(k))-LZheat)/csub.dzc(k);
        q[i][j][k] = heatfunc(csub.xc(i),csub.yc(j),
                              1.0-csub[i][j][k],csub.dV(i,j,k));
      } else {
        q[i][j][k] = heatfunc(csub.xc(i),csub.yc(j),
                              1.0,csub.dV(i,j,k));
      }
    }
  }

  if(boil::load_backup("time.txt",ts,time,
                       load_scalars, load_scalar_names,
                       load_vectors, load_vector_names,
                       load_nucls,   load_nucl_names,
                       load_cipcsls, load_cipcsl_names,
                       store_reals,  store_ints)) {
    conc.init();
    conc.color_to_vf();
  } else {

    /* start from single phase scratch */
    if(case_flag==0) {
#include "scratch_singlephase.cpp"
    /* create a bubble from initial time step */
    } else if(case_flag<5) {
#include "scratch_multiphase.cpp"
    } else {
      OMS(Underdevelopment. Exiting.);
      exit(0);
    }

  }

  cht.init();

  /* iint */
  //int iint = time.current_time() / t_per_plot;
  boil::oout<<"iint= "<<iint<<"\n";
