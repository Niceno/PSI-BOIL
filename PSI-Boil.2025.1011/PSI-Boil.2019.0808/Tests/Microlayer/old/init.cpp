  int ts;
  /* load variables */
  std::vector<Scalar*> load_scalars = { &press, &c.coarse };
  if(case_flag>2)
    load_scalars.push_back(&tpr.fine);
  else
    load_scalars.push_back(&tpr.coarse);
  std::vector<std::string> load_scalar_names = { "press", "c", "tpr" };

  std::vector<Vector*> load_vectors = { &uvw.coarse };
  std::vector<std::string> load_vector_names = { "uvw" };

  const real heater_extent = 1.5e-3;
  const real outer_factor = 1.;//30e-2;//18.5e-2;

  auto heatfunc = [&](const real x, const real y, const real c, const real v) {
    real h = qsrc*c*v;
    real ext0 = 4.3e-3;
    real ext1 = 4.3e-3;
    real ext_tilde = 0.8e-3;
    real b0 = 1.;
    real a0 = -b0/ext0;
    real q_tilde = a0*ext_tilde+b0;
    real a1 = q_tilde/(ext_tilde-ext1);
    real b1 = -a1*ext1;

    h *= std::max(0.,std::min(a0*x+b0,a1*x+b1));
  
    //h *= std::max(0.,(4.0e-3-x)/4.0e-3);
    //h *= std::min(std::max(0.,a0*x+b0,std::max(0.,a*x+b));
    //if(fabs(x)<=heater_extent&&fabs(y)<=heater_extent) {
    //} else {
    //  h *= outer_factor;
    //}
    return h;
  };

  /* solid */
  csub.fine   = 0.0;
  csub.coarse = 0.0;
  for_coarsefine(l) {
    for_vijk(csub[l],i,j,k) {
      if(csub[l].zc(k)<0.) {
        if(csub[l].zn(k+1)<=-LZheat) {
          csub[l][i][j][k] = 1.0;
        } else if(csub[l].zn(k)<= -LZheat) {
          csub[l][i][j][k] = (fabs(csub[l].zn(k))-LZheat)/csub[l].dzc(k);
          q[l][i][j][k] = heatfunc(csub[l].xc(i),csub[l].yc(j),
                                   1.0-csub[l][i][j][k],csub[l].dV(i,j,k));
        } else {
          q[l][i][j][k] = heatfunc(csub[l].xc(i),csub[l].yc(j),
                                   1.0,csub[l].dV(i,j,k));
        }
      }
    }
  }

  if(boil::load_backup("time.txt",ts,time,
                       load_scalars, load_scalar_names,
                       load_vectors, load_vector_names)) {
    conc_coarse.init();
    boil::prolongate_color_XZ(conc_coarse,conc_fine);
    conc_fine.color_to_vf();
    conc_fine.init();
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

  cht_fine.init();
  cht_coarse.init();

  /* set iint */
  int iint = time.current_time() / t_per_plot;
  boil::oout<<"iint= "<<iint<<"\n";
