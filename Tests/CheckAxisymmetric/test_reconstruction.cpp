static const real blending_angle = 40./180.*boil::pi;
static const real n0square = 1.-1./(1.+tan(blending_angle)*tan(blending_angle));
static const real n0 = sqrt(n0square);

#if 0
/*****************************************************************************/
void test_reconstruction_circle_xz(VOF & conc, Scalar & c, Scalar & ctest,
                                   Scalar & kappa, bool inverted, real radius,
                                   real xcent, real zcent,
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
      real theta = atan(fabs((ctest.zc(k)-zcent)/(ctest.xc(i)-xcent+boil::atto)));
      real nxtest = cos(theta);
      real nztest = sin(theta);

      if(ctest.xc(i)<xcent) {
        nxtest = -nxtest;
      }
      if(ctest.zc(k)<zcent) {
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
      real err = fabs(kap-kappa_real)/fabs(kappa_real);
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
#endif

/*****************************************************************************/
void test_reconstruction_sphere(bool partial,
                                VOFaxisym & conc, Scalar & c, Scalar & ctest,
                                Scalar & kappa, bool inverted, real radius,
                                real xcent, real ycent, real zcent,
                                NormMethod & nm,
                                std::vector<real> & nl1,
                                std::vector<real> & nli,
                                std::vector<real> & kl1,
                                std::vector<real> & kli) {

  conc.forward_cartesian(ctest);

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
      real x = fabs(ctest.xc(i)-xcent);
      real y = fabs(ctest.yc(j)-ycent);
      real z = fabs(ctest.zc(k)-zcent);
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

      if(ctest.xc(i)<xcent) {
        nxtest = -nxtest;
      }
      if(ctest.yc(j)<ycent) {
        nytest = -nytest;
      }
      if(ctest.zc(k)<zcent) {
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
      if(partial) {
        kappa_real /= 2.;
      }
      if(inverted)
        kappa_real = -kappa_real;
      real err = fabs(kap-kappa_real)/fabs(kappa_real);
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
void test_reconstruction_cone(const int NX, const int NZ,
                              VOFaxisym & conc, Scalar & c, Scalar & ctest,
                              Scalar & kappa, bool inverted,
                              real mx, real mz, real nalp,
                              NormMethod & nm,
                              std::vector<real> & nl1,
                              std::vector<real> & nli,
                              std::vector<real> & kl1,
                              std::vector<real> & kli) {

  conc.forward_cartesian(ctest);

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
    if(boil::realistic(kap)&&i>boil::BW+3&&k>boil::BW+3&&i<(NX+boil::BW-1-3)&&k<(NZ+boil::BW-1-3)) {
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
#if 0
      if(fabs(nxtest)>fabs(nztest)) {
        x0 = (nalp-fabs(nztest*ctest.zc(k)))/fabs(nxtest);
      } else {
        x0 = ctest.xc(i);
      }
#else 
      if       (fabs(nxtest)<n0) {
        x0 = ctest.xc(i);
      } else if(fabs(nztest)<n0) {
        x0 = (nalp-fabs(nztest*ctest.zc(k)))/fabs(nxtest);
      } else {
        real bfactor = (nztest*nztest - n0square)/(1.-2.*n0square);
  #if 0
        real bangle = atan(sqrt(1./(nztest*nztest)-1.));
        boil::oout<<nztest<<" "<<nztest*nztest-n0square<<" "<<(1.-2.*n0square)<<" | "<<bfactor<<" "<<(1./(1.+pow(tan(bangle),2.))-n0square)/(1.-2.*n0square)<<" "<<nztest<<" "<<bangle/boil::pi*180.<<boil::endl;
  #endif
        x0 = ctest.xc(i)*bfactor
           + (nalp-fabs(nztest*ctest.zc(k)))/fabs(nxtest) * (1.-bfactor);
      }
#endif
      real kappa_real = nxtest/x0;
      real err = fabs((kap-kappa_real)/kappa_real);
      if(err>err_k_linf) {
        err_k_linf = err;
      }
      err_k_l1 += err;

      boil::oout<<i<<" "<<k<<" | "<<" "<<x0<<" "<<kap<<" "<<kappa_real<<" | "<<err<<" "<<err_k_l1<<" "<<err_k_linf<<boil::endl;
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
void test_reconstruction_doughnut(bool partial_major, bool partial_minor,
                                 VOFaxisym & conc, Scalar & c, Scalar & ctest,
                                 Scalar & kappa, bool inverted, 
                                 real radius_major, real radius_minor,
                                 real xcent, real zcent,
                                 NormMethod & nm,
                                 std::vector<real> & nl1,
                                 std::vector<real> & nli,
                                 std::vector<real> & kl1,
                                 std::vector<real> & kli) {

  conc.forward_cartesian(ctest);

  real err_n_l1(0.0),  err_n_linf(0.0);
  real err_c_l1(0.0),  err_c_linf(0.0);
  int cnt(0);
  for_vijk(ctest,i,j,k) {
    if(conc.adens[i][j][k]>boil::pico) {
      cnt++;
      real theta = atan(fabs((ctest.zc(k)-zcent)/(ctest.xc(i)-xcent+boil::atto)));
      real nxtest = cos(theta);
      real nztest = sin(theta);

      if(ctest.xc(i)<xcent) {
        nxtest = -nxtest;
      }
      if(ctest.zc(k)<zcent) {
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
    //if(boil::realistic(kap)&&i>boil::BW+3&&i<(NX+boil::BW-1-3)) {
      knt++;
      real kappa_major(0.0);
      real kappa_minor(0.0);
      if(partial_minor) kappa_minor = 1./radius_minor;
      if(partial_major) {

        real theta = atan(fabs((ctest.zc(k)-zcent)/(ctest.xc(i)-xcent+boil::atto)));
        real nxtest = fabs(cos(theta));
        real nztest = fabs(sin(theta));

        real x0;
        if       (fabs(nxtest)<n0) {
          x0 = ctest.xc(i)-xcent;
        } else if(fabs(nztest)<n0) {
          real z0 = ctest.zc(k)-zcent;
          x0 = sqrt(radius_minor*radius_minor-z0*z0);
        } else {
          real z0 = ctest.zc(k)-zcent;
          real bfactor = (nztest*nztest - n0square)/(1.-2.*n0square);
          x0 = (ctest.xc(i)-zcent)*bfactor
             + sqrt(radius_minor*radius_minor-z0*z0) * (1.-bfactor);
        }


        real xpos = (ctest.xc(i)-xcent);
        kappa_major = xpos/radius_minor/(radius_major+xpos);
      }
      if(inverted) {
        kappa_minor = -kappa_minor;
        kappa_major = -kappa_major;
      }
      real kappa_real = kappa_major + kappa_minor;
      real err = fabs(kap-kappa_real)/fabs(kappa_real);
      //real err = fabs(kap-kappa_real);
      if(err>err_k_linf) {
        err_k_linf = err;
      }
      err_k_l1 += err;
      
      boil::oout<<i<<" "<<k<<" | "<<" "<<radius_major<<" "<<radius_minor<<" "<<kap<<" "<<kappa_real<<" | "<<err<<" "<<err_k_l1<<" "<<err_k_linf<<boil::endl;
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
