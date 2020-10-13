  int ts;
  /* load variables */
  std::vector<Scalar*> load_scalars = { &press, &c.coarse };
  if(case_flag==2)
    load_scalars.push_back(&tpr.fine);
  else
    load_scalars.push_back(&tpr.coarse);
  std::vector<std::string> load_scalar_names = { "press", "c", "tpr" };

  std::vector<Vector*> load_vectors = { &uvw.coarse };
  std::vector<std::string> load_vector_names = { "uvw" };

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
          /* estimation */ 
          if(fabs(csub[l].xc(i))<=0.3*LX1&&fabs(csub[l].yc(j))<=0.3*LX1)
            q[l][i][j][k] = qsrc*(1.0-csub[l][i][j][k])*csub[l].dV(i,j,k);
          else
            q[l][i][j][k] = 0.6*qsrc*(1.0-csub[l][i][j][k])*csub[l].dV(i,j,k);
        } else {
          if(fabs(csub[l].xc(i))<=0.3*LX1&&fabs(csub[l].yc(j))<=0.3*LX1)
            q[l][i][j][k] = qsrc*csub[l].dV(i,j,k);
          else
            q[l][i][j][k] = 0.6*qsrc*csub[l].dV(i,j,k);
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

  /* set iint */
  int iint = time.current_time() / t_per_plot;
  boil::oout<<"iint= "<<iint<<"\n";
