/*****************************************************************************/
void test_reconstruction_circle_xz(VOF & conc, Scalar & c, Scalar & ctest,
                                   Scalar & kappa, bool inverted, real radius,
                                   NormMethod & nm,
                                   std::vector<real> & nl1,
                                   std::vector<real> & nli,
                                   std::vector<real> & kl1,
                                   std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l1(0.0),  err_n_linf(0.0);
  real err_c_l1(0.0),  err_c_linf(0.0);
  int cnt(0);
  for_vijk(ctest,i,j,k) {
    if(conc.adens[i][j][k]>boil::pico) {
      cnt++;
      real theta = atan(fabs(ctest.zc(k)/(ctest.xc(i)+boil::atto)));
      real nxtest = cos(theta);
      real nztest = sin(theta);

      if(ctest.xc(i)<0) {
        nxtest = -nxtest;
      }
      if(ctest.zc(k)<0) {
        nztest = -nztest;
      }
      if(!inverted) {
        nxtest = -nxtest;
        nztest = -nztest;
      }

      real nnx = fabs( (conc.nx)[i][j][k] - nxtest);
      real nnz = fabs( (conc.nz)[i][j][k] - nztest);

      real err = sqrt(nnx*nnx + nnz*nnz);
      if(err>err_n_linf) {
        err_n_linf = err;
      }
      err_n_l1 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l1 += cc;
    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l1);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l1);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l1 /= real(cnt);
    err_c_l1 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l1<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l1<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l1(0.0),  err_k_linf(0.0);
  int knt(0);
  for_vijk(ctest,i,j,k) {
    real kap = kappa[i][j][k];
    if(boil::realistic(kap)) {
      knt++;
      real kappa_real = 1./radius;
      if(inverted)
        kappa_real = -kappa_real;
      real err = fabs(kap-kappa_real);
      if(err>err_k_linf) {
        err_k_linf = err;
      }
      err_k_l1 += err;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l1);

  if(knt>0) {
    err_k_l1 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l1<<" (in "<<knt<<" cells)"<<boil::endl;

  nl1.push_back(err_n_l1);
  nli.push_back(err_n_linf);
  kl1.push_back(err_k_l1);
  kli.push_back(err_k_linf);

  return;
}

/*****************************************************************************/
void test_reconstruction_triangle_xz(VOF & conc, Scalar & c, Scalar & ctest,
                                     Scalar & kappa, bool inverted,
                                     real mx, real mz,
                                     NormMethod & nm,
                                     std::vector<real> & nl1,
                                     std::vector<real> & nli,
                                     std::vector<real> & kl1,
                                     std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l1(0.0),  err_n_linf(0.0);
  real err_c_l1(0.0),  err_c_linf(0.0);
  int cnt(0);
  for_vijk(ctest,i,j,k) {
    if(conc.adens[i][j][k]>boil::pico&&i>boil::BW&&k>boil::BW) {
      cnt++;
      real nxtest = mx;
      real nztest = mz;

      if(ctest.xc(i)<0) {
        nxtest = -nxtest;
      }
      if(ctest.zc(k)<0) {
        nztest = -nztest;
      }
      if(!inverted) {
        nxtest = -nxtest;
        nztest = -nztest;
      }

      real nnx = fabs( (conc.nx)[i][j][k] - nxtest);
      real nnz = fabs( (conc.nz)[i][j][k] - nztest);

      real err = sqrt(nnx*nnx + nnz*nnz);
      if(err>err_n_linf) {
        err_n_linf = err;
      }
      err_n_l1 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l1 += cc;

    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l1);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l1);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l1 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l1<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l1<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l1(0.0),  err_k_linf(0.0);
  int knt(0);
  for_vijk(ctest,i,j,k) {
    real kap = kappa[i][j][k];
    if(boil::realistic(kap)&&i>boil::BW+3&&k>boil::BW+3&&conc.adens[i][j][k]>boil::pico) {
      knt++;
      real kappa_real = 0.0;
      if(inverted)
        kappa_real = -kappa_real;
      real err = fabs(kap-kappa_real);
      if(err>err_k_linf) {
        err_k_linf = err;
      }
      err_k_l1 += err;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l1);

  if(knt>0) {
    err_k_l1 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l1<<" (in "<<knt<<" cells)"<<boil::endl;

  nl1.push_back(err_n_l1);
  nli.push_back(err_n_linf);
  kl1.push_back(err_k_l1);
  kli.push_back(err_k_linf);

  return;
}

/*****************************************************************************/
void test_reconstruction_sphere(VOF & conc, Scalar & c, Scalar & ctest,
                                Scalar & kappa, bool inverted, real radius,
                                NormMethod & nm,
                                std::vector<real> & nl1,
                                std::vector<real> & nli,
                                std::vector<real> & kl1,
                                std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l1(0.0),  err_n_linf(0.0);
  real err_c_l1(0.0),  err_c_linf(0.0);
  int cnt(0);
  for_vijk(ctest,i,j,k) {
    if(conc.adens[i][j][k]>boil::pico) {
      cnt++;

      /* surface normal to a sphere
         theta = inclination, phi = azimuth
         x = r*sin(theta)*cos(phi)
         y = r*sin(theta)*sin(phi)
         z = r*cos(theta) */
      real x = fabs(ctest.xc(i));
      real y = fabs(ctest.yc(j));
      real z = fabs(ctest.zc(k));
      real r = sqrt(x*x+y*y+z*z);
      real phi = atan(y/(x+boil::atto));
      real theta = acos(z/(r+boil::atto));

      real nxtest = sin(theta)*cos(phi);
      real nytest = sin(theta)*sin(phi);
      real nztest = cos(theta);

      real nsum = sqrt(nxtest*nxtest+nytest*nytest+nztest*nztest+boil::atto);
      nxtest /= nsum;
      nytest /= nsum;
      nztest /= nsum;

      if(ctest.xc(i)<0) {
        nxtest = -nxtest;
      }
      if(ctest.yc(j)<0) {
        nytest = -nytest;
      }
      if(ctest.zc(k)<0) {
        nztest = -nztest;
      }
      if(!inverted) {
        nxtest = -nxtest;
        nytest = -nytest;
        nztest = -nztest;
      }

      real nnx = fabs( (conc.nx)[i][j][k] - nxtest);
      real nny = fabs( (conc.ny)[i][j][k] - nytest);
      real nnz = fabs( (conc.nz)[i][j][k] - nztest);

      real err = sqrt(nnx*nnx + nny*nny + nnz*nnz);
      if(err>err_n_linf) {
        err_n_linf = err;
        //boil::oout<<i<<" "<<j<<" "<<k<<" | "<<(conc.nx)[i][j][k]<<" "<<(conc.ny)[i][j][k]<<" "<<(conc.nz)[i][j][k]<<" | "<<nxtest<<" "<<nytest<<" "<<nztest<<boil::endl;
      }
      err_n_l1 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l1 += cc;
    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l1);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l1);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l1 /= real(cnt);
    err_c_l1 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l1<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l1<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l1(0.0),  err_k_linf(0.0);
  int knt(0);
  for_vijk(ctest,i,j,k) {
    real kap = kappa[i][j][k];
    if(boil::realistic(kap)) {
      knt++;
      real kappa_real = 2./radius;
      if(inverted)
        kappa_real = -kappa_real;
      real err = fabs(kap-kappa_real);
      if(err>err_k_linf) {
        err_k_linf = err;
      }
      err_k_l1 += err;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l1);
  boil::cart.sum_int(&knt);

  if(knt>0) {
    err_k_l1 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l1<<" (in "<<knt<<" cells)"<<boil::endl;

  nl1.push_back(err_n_l1);
  nli.push_back(err_n_linf);
  kl1.push_back(err_k_l1);
  kli.push_back(err_k_linf);

  return;
}

/*****************************************************************************/
void test_reconstruction_cone(VOF & conc, Scalar & c, Scalar & ctest,
                              Scalar & kappa, bool inverted,
                              real mx, real mz, real nalp,
                              NormMethod & nm,
                              std::vector<real> & nl1,
                              std::vector<real> & nli,
                              std::vector<real> & kl1,
                              std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l1(0.0),  err_n_linf(0.0);
  real err_c_l1(0.0),  err_c_linf(0.0);
  int cnt(0);
  for_vijk(ctest,i,j,k) {
    if(conc.adens[i][j][k]>boil::pico&&i>boil::BW&&k>boil::BW) {
      cnt++;
      real nxtest = mx;
      real nztest = mz;

      if(ctest.xc(i)<0) {
        nxtest = -nxtest;
      }
      if(ctest.zc(k)<0) {
        nztest = -nztest;
      }
      if(!inverted) {
        nxtest = -nxtest;
        nztest = -nztest;
      }

      real nnx = fabs( (conc.nx)[i][j][k] - nxtest);
      real nnz = fabs( (conc.nz)[i][j][k] - nztest);

      real err = sqrt(nnx*nnx + nnz*nnz);
      if(err>err_n_linf) {
        err_n_linf = err;
      }
      err_n_l1 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l1 += cc;

    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l1);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l1);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l1 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l1<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l1<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l1(0.0),  err_k_linf(0.0);
  int knt(0);
  for_vijk(ctest,i,j,k) {
    real kap = kappa[i][j][k];
    if(boil::realistic(kap)&&i>boil::BW+3&&k>boil::BW+3) {
      knt++;
      real nxtest = mx;
      real nztest = mz;
      if(ctest.zc(k)<0) {
        nztest = -nztest;
      }
      if(ctest.xc(i)<0) {
        nxtest = -nxtest;
      }
      if(inverted) {
        nxtest = -nxtest;
        nztest = -nztest;
      }
      real x0;
      if(fabs(nxtest)>fabs(nztest)) {
        x0 = (nalp-fabs(nztest*ctest.zc(k)))/fabs(nxtest);
      } else {
        x0 = ctest.xc(i);
      }
      real kappa_real = nxtest/x0;
      real err = fabs((kap-kappa_real)/kappa_real);
      if(err>err_k_linf) {
        err_k_linf = err;
      }
      err_k_l1 += err;

      boil::oout<<i<<" "<<k<<" "<<kap<<" "<<kappa_real<<" | "<<err<<" "<<err_k_l1<<" "<<err_k_linf<<boil::endl;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l1);

  if(knt>0) {
    err_k_l1 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l1<<" (in "<<knt<<" cells)"<<boil::endl;

  nl1.push_back(err_n_l1);
  nli.push_back(err_n_linf);
  kl1.push_back(err_k_l1);
  kli.push_back(err_k_linf);

  return;
}
