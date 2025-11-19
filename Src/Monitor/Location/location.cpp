#include "location.h"

/***************************************************************************//**
*  \brief Constructor for monitoring Location in computational domain
*
*  \param n - Monitor's name (it is printed before the value),
*  \param d - Domain on which the monitor is created,
*  \param i, j, k - logical coordinates defining Location's position.
*
*  If location can not be found in prescribed domain, it writes the error
*  message and stops the program.
*******************************************************************************/
Location::Location(const char * n, 
                   const Domain & d, const int I, const int J, const int K) 
 : name(n), dom(&d), m_i(I+boil::BW-1), m_j(J+boil::BW-1), m_k(K+boil::BW-1) {

  assert(true ==1);
  assert(false==0);

  found_here = dom->contains_IJK(m_i, m_j, m_k); // true=1; false=0

  int found_anywhere = found_here;
  boil::cart.sum_int(&found_anywhere);

  /* if monitoring point not found, exit */
  if( !found_anywhere ) {
    boil::oout << "Location " << name << " couldn't be found " << boil::endl;
    exit(0);
  }

  /* monitoring point found, print it */
  if( found_here ) {
    dom->locals(&m_i, &m_j, &m_k); // send global, return local
    if(name)
    boil::aout << "Location " << name <<
                  " at x = " << dom->xc(m_i) <<
                    ", y = " << dom->yc(m_j) <<
                    ", z = " << dom->zc(m_k) << " created." << boil::endl;
  }
}

/******************************************************************************/
void Location::print(const Scalar & phi) {

  real value=0;
  real count=0;
  real x=0, y=0, z=0;

  if( found_here ) {
    value = phi[m_i][m_j][m_k];
    count = 1;
    x=phi.xc(m_i);
    y=phi.yc(m_j);
    z=phi.zc(m_k);
  }

  boil::cart.sum_real(&value);
  boil::cart.sum_real(&count);
  boil::cart.sum_real(&x);
  boil::cart.sum_real(&y);
  boil::cart.sum_real(&z);

  assert(count > 0.0);

  if(name) {
    boil::oout << name << ": " << value/count << boil::endl;
  } else {
    boil::oout << x/count << " " <<
                  y/count << " " <<
                  z/count << " " << value/count << boil::endl;
  }
}

/******************************************************************************/
void Location::print(const Vector & u, const Comp & m) {

  real value=0;
  real count=0;
  real x=0, y=0, z=0;

  if( found_here ) {
    value = u[m][m_i][m_j][m_k];
    count = 1;
    x=u.xc(m,m_i);
    y=u.yc(m,m_j);
    z=u.zc(m,m_k);
  }

  boil::cart.sum_real(&value);
  boil::cart.sum_real(&count);
  boil::cart.sum_real(&x);
  boil::cart.sum_real(&y);
  boil::cart.sum_real(&z);

  assert(count > 0.0);

  if(name) {
    boil::oout << name << ": " << value/count << boil::endl;
  } else {
    boil::oout << x/count << " " <<
                  y/count << " " <<
                  z/count << " " << value/count << boil::endl;
  }
}

/******************************************************************************/
real Location::value(const Scalar & phi) {

  real val=0;
  real count=0;

  if( found_here ) {
    val = phi[m_i][m_j][m_k];
    count = 1.0;
  }

  boil::cart.sum_real(&val);
  boil::cart.sum_real(&count);

  assert(count > 0.0);

  return val/count;
}

/******************************************************************************/
real Location::value(const Vector & u, const Comp & m) {

  real val=0;
  real count=0;

  if( found_here ) {
    val = u[m][m_i][m_j][m_k];
    count = 1;
  }

  boil::cart.sum_real(&val);
  boil::cart.sum_real(&count);

  assert(count > 0.0);

  return val/count;
}

/******************************************************************************/
real Location::get_scalar(const Scalar & phi) {

  real value=0;
  real count=0;

  if( found_here ) {
    value = phi[m_i][m_j][m_k];
    count = 1;
  }

  boil::cart.sum_real(&value);
  boil::cart.sum_real(&count);

  assert(count > 0.0);

  return value/count;
}

/******************************************************************************/
real Location::get_vector(const Vector & u, const Comp & m) {

  real value=0;
  real count=0;

  if( found_here ) {
    value = u[m][m_i][m_j][m_k];
    count = 1;
  }

  boil::cart.sum_real(&value);
  boil::cart.sum_real(&count);

  assert(count > 0.0);

  return value/count;
}

/******************************************************************************/
real Location::get_grid(const Scalar & phi, const Comp & m) {

  real count=0;
  real xyz=0;

  if( found_here ) {
    count = 1;
    if (m==Comp::x()) {
      xyz=phi.xc(m_i);
    } else if (m==Comp::y()) {
      xyz=phi.yc(m_j);
    } else {
      xyz=phi.zc(m_k);
    }
  }

  boil::cart.sum_real(&count);
  boil::cart.sum_real(&xyz);

  assert(count > 0.0);

  return xyz/count;
}

/******************************************************************************/
real Location::get_grid(const Vector & u, const Comp & m, const Comp & n) {

  real count=0;
  real xyz=0;

  if( found_here ) {
    count = 1;
    if (n==Comp::x()) {
      xyz=u.xc(m,m_i);
    } else if (n==Comp::y()) {
      xyz=u.yc(m,m_j);
    } else {
      xyz=u.zc(m,m_k);
    }
  }

  boil::cart.sum_real(&count);
  boil::cart.sum_real(&xyz);

  assert(count > 0.0);

  return xyz/count;
}

