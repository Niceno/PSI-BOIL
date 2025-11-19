#include "scalar.h"

/******************************************************************************/
real Scalar::min() const {
  real vmin = val[boil::BW][boil::BW][boil::BW];
  for_ijk(i,j,k) {
    if( val[i][j][k] < vmin ) {
      vmin = val[i][j][k];
    }
  }
  boil::cart.min_real(&vmin);
  return vmin;
}

/******************************************************************************/
real Scalar::min_abs() const {
  real vmin = fabs(val[boil::BW][boil::BW][boil::BW]);
  for_ijk(i,j,k) {
    if( fabs(val[i][j][k]) < vmin ) {
      vmin = fabs(val[i][j][k]);
    }
  }
  boil::cart.min_real(&vmin);
  return vmin;
}

/******************************************************************************/
real Scalar::min_at() const {
  real vmin = val[boil::BW][boil::BW][boil::BW];
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k] < vmin ) {
      vmin = val[i][j][k]; im=i; jm=j; km=k;
    }
  }
  real vmin_l = vmin;
  boil::cart.min_real(&vmin);
  if(vmin == vmin_l) {
    boil::aout << name() << " min = " << vmin << " at: " 
               << xc(im) << " " 
               << yc(jm) << " " 
               << zc(km) << boil::endl;
  }
  return vmin;
}

/******************************************************************************/
real Scalar::max() const {
  real vmax = val[boil::BW][boil::BW][boil::BW];
  for_ijk(i,j,k) {
    if( val[i][j][k] > vmax ) {
      vmax = val[i][j][k];
    }
  }
  boil::cart.max_real(&vmax);
  return vmax;
}

/******************************************************************************/
real Scalar::max_abs() const {
  real vmax = fabs(val[boil::BW][boil::BW][boil::BW]);
  for_ijk(i,j,k) {
    if( fabs(val[i][j][k]) > vmax ) {
      vmax = fabs(val[i][j][k]);
    }
  }
  boil::cart.max_real(&vmax);
  return vmax;
}

/******************************************************************************/
real Scalar::max_at() const {
  real vmax = val[boil::BW][boil::BW][boil::BW];
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k] > vmax ) {
      vmax = val[i][j][k]; im=i; jm=j; km=k;
    }
  }
  real vmax_l = vmax;
  boil::cart.max_real(&vmax);
  if(vmax == vmax_l) {
    boil::aout << name() << " max = " << vmax << " at: " 
               << xc(im) << " " 
               << yc(jm) << " " 
               << zc(km) << boil::endl;
  }
  return vmax;
}

/******************************************************************************/
real Scalar::min_voldiv() const {
  real vmin = val[boil::BW][boil::BW][boil::BW]/dV(boil::BW,boil::BW,boil::BW);
  for_ijk(i,j,k) {
    if( val[i][j][k]/dV(i,j,k) < vmin ) {
      vmin = val[i][j][k]/dV(i,j,k);
    }
  }
  boil::cart.min_real(&vmin);
  return vmin;
}

/******************************************************************************/
real Scalar::min_abs_voldiv() const {
  real vmin = fabs(val[boil::BW][boil::BW][boil::BW]
                   /dV(boil::BW,boil::BW,boil::BW));
  for_ijk(i,j,k) {
    if( fabs(val[i][j][k]/dV(i,j,k)) < vmin ) {
      vmin = fabs(val[i][j][k]/dV(i,j,k));
    }
  }
  boil::cart.min_real(&vmin);
  return vmin;
}

/******************************************************************************/
real Scalar::min_at_voldiv() const {
  real vmin = val[boil::BW][boil::BW][boil::BW]/dV(boil::BW,boil::BW,boil::BW);
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k]/dV(i,j,k) < vmin ) {
      vmin = val[i][j][k]/dV(i,j,k); im=i; jm=j; km=k;
    }
  }
  real vmin_l = vmin;
  boil::cart.min_real(&vmin);
  if(vmin == vmin_l) {
    boil::aout << name() << " min_voldiv = " << vmin << " at: "
               << xc(im) << " "
               << yc(jm) << " "
               << zc(km) << boil::endl;
  }
  return vmin;
}

/******************************************************************************/
real Scalar::max_voldiv() const {
  real vmax = val[boil::BW][boil::BW][boil::BW]/dV(boil::BW,boil::BW,boil::BW);
  for_ijk(i,j,k) {
    if( val[i][j][k]/dV(i,j,k) > vmax ) {
      vmax = val[i][j][k]/dV(i,j,k);
    }
  }
  boil::cart.max_real(&vmax);
  return vmax;
}

/******************************************************************************/
real Scalar::max_abs_voldiv() const {
  real vmax = fabs(val[boil::BW][boil::BW][boil::BW]
                   /dV(boil::BW,boil::BW,boil::BW));
  for_ijk(i,j,k) {
    if( fabs(val[i][j][k]/dV(i,j,k)) > vmax ) {
      vmax = fabs(val[i][j][k]/dV(i,j,k));
    }
  }
  boil::cart.max_real(&vmax);
  return vmax;
}

/******************************************************************************/
real Scalar::max_at_voldiv() const {
  real vmax = val[boil::BW][boil::BW][boil::BW]/dV(boil::BW,boil::BW,boil::BW);
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k]/dV(i,j,k) > vmax ) {
      vmax = val[i][j][k]/dV(i,j,k); im=i; jm=j; km=k;
    }
  }
  real vmax_l = vmax;
  boil::cart.max_real(&vmax);
  if(vmax == vmax_l) {
    boil::aout << name() << " max_voldiv = " << vmax << " at: "
               << xc(im) << " "
               << yc(jm) << " "
               << zc(km) << boil::endl;
  }
  return vmax;
}

