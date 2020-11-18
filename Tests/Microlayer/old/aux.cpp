/* reading settings from the input file */
template <class T>
bool read_variable(std::ifstream & input, T & var) {
  std::string line;
  if(getline(input,line)) {
    std::istringstream ssline(line);
    ssline >> var;
    OPR(var);
    return 0;
  } else {
    return 1;
  }
}

/* write an interpolation table from a single-proc scalar;
 * the vector is for solid-liquid boundary
 * 1st line = x-vals
 * 2nd line = z-vals
 * rest of the file data in the format val(z = col,x = row)
 */
void write_table(const Scalar & sca, const Vector & bndsca,
                 std::fstream & otp) {
  otp << 0. << " ";
  for_vi(sca,i) {
    otp << sca.xc(i) << " ";
  }
  otp << boil::endl;

  bool zerocrossed(0);
  for_vk(sca,k) {
    if(sca.zc(k)<0.) {
      otp << sca.zc(k) << " ";
    } else if(!zerocrossed) {
      otp << 0. << " " << sca.zc(k) << " ";
      zerocrossed = true;
    } else {
      otp << sca.zc(k) << " ";
    }
  }
  otp << boil::endl;

  int j(boil::BW);
  zerocrossed = false;
  for_vk(sca,k) {
    if(sca.zc(k)<0.) {
      otp<<sca[boil::BW][j][k]<<" ";
      for_vi(sca,i) {
        otp<<sca[i][j][k]<<" ";
      }
      otp << boil::endl;
    } else if(!zerocrossed) {
      otp<<bndsca[Comp::k()][boil::BW][j][k]<<" ";
      for_vi(sca,i) {
        otp<<bndsca[Comp::k()][i][j][k]<<" ";
      }
      otp << boil::endl;
      otp<<sca[boil::BW][j][k]<<" ";
      for_vi(sca,i) {
        otp<<sca[i][j][k]<<" ";
      }
      otp << boil::endl;
      zerocrossed = true;
    } else {
      otp<<sca[boil::BW][j][k]<<" ";
      for_vi(sca,i) {
        otp<<sca[i][j][k]<<" ";
      }
      otp << boil::endl;
    }
  }

  return;
}

/* output to a file */
void output_to_file(const Scalar & sca, const Vector & bndsca) {

  std::string fname = "temperature.input";
  std::fstream output;
  output.open(fname,std::ios::out);
  write_table(sca,bndsca,output);
  output.close();

  return;
}

/* read an interpolation table:
 * 1st line = x-vals
 * 2nd line = z-vals
 * rest of the file data in the format val(z = col,x = row)
 */
void read_table(std::vector< std::vector<real> > & table,
                std::fstream & inp) {

  assert(table.empty());

  /* buffer */
  std::string line;

  /* read x */
  std::getline(inp,line);
  std::istringstream issx(line);
  std::vector<real> xvals { std::istream_iterator<real>(issx),
                            std::istream_iterator<real>() };
  const int xsize = xvals.size();
  table.push_back(xvals);

  /* read z */
  std::getline(inp,line);
  std::istringstream issz(line);
  std::vector<real> zvals { std::istream_iterator<real>(issz),
                            std::istream_iterator<real>() };
  const int zsize = zvals.size();
  table.push_back(zvals);

  /* read the rest */
  while(std::getline(inp,line)) {
    std::istringstream iss(line);
    std::vector<real> vals { std::istream_iterator<real>(iss),
                             std::istream_iterator<real>() };
    assert(vals.size() == xsize);
    table.push_back(vals);
  }
  assert( (table.size()-2) == zsize );

  return;
}

/* bilinear interpolation from a table */
real bilinear(const real xx, const real zz, 
              const std::vector< std::vector<real> > & table) {
  real x(xx),z(zz);
  int l0(-1),m0(-1);
  real x0(-1.),x1(-1.);
  real z0(-1.),z1(-1.);
  if(x<table[0][0]) {
    x = table[0][0];
  } else if(x>table[0][table[0].size()-1]) {
    x = table[0][table[0].size()-1];
  }
  for(int l(0); l<table[0].size()-1; ++l) {
    if(table[0][l]<=x && x<=table[0][l+1]) {
      l0 = l;
      x0 = table[0][l];
      x1 = table[0][l+1];
      break;
    }
  }
  if(z<table[1][0]) {
    z = table[1][0];
  } else if(z>table[1][table[1].size()-1]) {
    z = table[1][table[1].size()-1];
  }
  for(int m(0); m<table[1].size()-1; ++m) {
    if(table[1][m]<=z && z<=table[1][m+1]) {
      m0 = m+2; /* shift due to the first rows being x,z */
      z0 = table[1][m];
      z1 = table[1][m+1];
      break;
    }
  }

  return (z1-z)/(z1-z0) * ( (x1-x)/(x1-x0) * table[m0  ][l0  ]
                           +(x-x0)/(x1-x0) * table[m0  ][l0+1])
        +(z-z0)/(z1-z0) * ( (x1-x)/(x1-x0) * table[m0+1][l0  ]
                           +(x-x0)/(x1-x0) * table[m0+1][l0+1]);
}

/* interpolation from a file */
bool interpolate_from_file(Scalar & sca) {

  /* reading of the interpolation table done sequentially
   * (it is possible to implement parallel reading though) */
  std::string fname = "temperature.input";
  std::fstream input;
  std::vector< std::vector<real> > table;
  int open_failed(0);

  for(int pid(0); pid<boil::cart.nproc(); ++pid) {
    if(pid==boil::cart.iam()) {
      input.open(fname,std::ios::in);
      open_failed = input.fail();
      if(!open_failed) {
        read_table(table,input);
      }
      input.close();
    }
    boil::cart.barrier();
  }
  boil::cart.sum_int(&open_failed);

  if(open_failed) {
    return 1;
  } else {
    /* bilinear interpolation */
    for_vijk(sca,i,j,k) {
      sca[i][j][k] = bilinear(sca.xc(i),sca.zc(k),table);
    }
  }

  return 0;
}

/* row output */
real output_row(const real zp, const Scalar & sca, const bool inv) {
  real h0 = 0.0;
  for_vk(sca,k) {
    if(sca.zn(k)<zp&&sca.zn(k+1)>=zp) {
      for_vi(sca,i) {
        h0 += (inv ? (1.-sca[i][boil::BW][k]) : sca[i][boil::BW][k]) * sca.dxc(i);
      } /* i */
    } /* in range */
  } /* k */
  boil::cart.sum_real(&h0);
  return h0;
}
