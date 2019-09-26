void test_reconstruction_circle_xz(VOF & conc, Scalar & c, Scalar & ctest,
                                   Scalar & kappa, bool inverted, real radius,
                                   NormMethod & nm,
                                   std::vector<real> & nl2,
                                   std::vector<real> & nli,
                                   std::vector<real> & kl2,
                                   std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l2(0.0),  err_n_linf(0.0);
  real err_c_l2(0.0),  err_c_linf(0.0);
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
      err_n_l2 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l2 += cc;
    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l2);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l2);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l2 /= real(cnt);
    err_c_l2 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l2<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l2<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l2(0.0),  err_k_linf(0.0);
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
      err_k_l2 += err;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l2);

  if(knt>0) {
    err_k_l2 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l2<<" (in "<<knt<<" cells)"<<boil::endl;

  nl2.push_back(err_n_l2);
  nli.push_back(err_n_linf);
  kl2.push_back(err_k_l2);
  kli.push_back(err_k_linf);

  return;
}

void test_reconstruction_triangle_xz(VOF & conc, Scalar & c, Scalar & ctest,
                                     Scalar & kappa, bool inverted,
                                     real mx, real mz,
                                     NormMethod & nm,
                                     std::vector<real> & nl2,
                                     std::vector<real> & nli,
                                     std::vector<real> & kl2,
                                     std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l2(0.0),  err_n_linf(0.0);
  real err_c_l2(0.0),  err_c_linf(0.0);
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
      err_n_l2 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l2 += cc;

    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l2);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l2);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l2 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l2<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l2<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l2(0.0),  err_k_linf(0.0);
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
      err_k_l2 += err;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l2);

  if(knt>0) {
    err_k_l2 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l2<<" (in "<<knt<<" cells)"<<boil::endl;

  nl2.push_back(err_n_l2);
  nli.push_back(err_n_linf);
  kl2.push_back(err_k_l2);
  kli.push_back(err_k_linf);

  return;
}

void test_reconstruction_sphere(VOF & conc, Scalar & c, Scalar & ctest,
                                Scalar & kappa, bool inverted, real radius,
                                NormMethod & nm,
                                std::vector<real> & nl2,
                                std::vector<real> & nli,
                                std::vector<real> & kl2,
                                std::vector<real> & kli) {

  conc.forward(ctest);

  real err_n_l2(0.0),  err_n_linf(0.0);
  real err_c_l2(0.0),  err_c_linf(0.0);
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
      err_n_l2 += err;

      real cc = fabs(ctest[i][j][k]-c[i][j][k]);
      if(cc>err_c_linf) {
        err_c_linf = cc;
      }
      err_c_l2 += cc;
    }
  }
  boil::cart.max_real(&err_n_linf);
  boil::cart.sum_real(&err_n_l2);
  boil::cart.max_real(&err_c_linf);
  boil::cart.sum_real(&err_c_l2);
  boil::cart.sum_int(&cnt);

  if(cnt>0) {
    err_n_l2 /= real(cnt);
    err_c_l2 /= real(cnt);
  }

  boil::oout<<"Error-normal-vector-for "<<nm<<": "<<err_n_linf<<" "<<err_n_l2<<" (in "<<cnt<<" cells)"<<boil::endl;
  boil::oout<<"Error-color-for "<<nm<<": "<<err_c_linf<<" "<<err_c_l2<<" (in "<<cnt<<" cells)"<<boil::endl;

  /* curvature */
  real err_k_l2(0.0),  err_k_linf(0.0);
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
      err_k_l2 += err;
    }
  }
  boil::cart.max_real(&err_k_linf);
  boil::cart.sum_real(&err_k_l2);
  boil::cart.sum_int(&knt);

  if(knt>0) {
    err_k_l2 /= real(knt);
  }

  boil::oout<<"Error-kappa-for "<<nm<<": "<<err_k_linf<<" "<<err_k_l2<<" (in "<<knt<<" cells)"<<boil::endl;

  nl2.push_back(err_n_l2);
  nli.push_back(err_n_linf);
  kl2.push_back(err_k_l2);
  kli.push_back(err_k_linf);

  return;
}

