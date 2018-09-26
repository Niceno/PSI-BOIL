#include "rack.h"

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n   - rack's name (it is printed before the values),
*  \param d   - Domain on which the rack is created,
*  \param ri  - range of i's (minimum and maximum)
*  \param j,k - logical coordinates defining racks's position.
*
*  \note This constructor creates the rack in "i" direction.                     
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d, 
           const Range<int> ri, const int j, const int k) 
 : name(n), r_i(ri), r_j(j,j), r_k(k,k) {

  int size = ri.last() - ri.first() + 1;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) {
    mons[c++] = new Location(NULL, d, i, j, k);
  }
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n   - rack's name (it is printed before the values),
*  \param d   - Domain on which the rack is created,
*  \param rj  - range of i's (minimum and maximum)
*  \param i,k - logical coordinates defining racks's position.
*
*  \note This constructor creates the rack in "j" direction.                     
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d, 
           const int i, const Range<int> rj, const int k) 
 : name(n), r_i(i,i), r_j(rj), r_k(k,k) {

  int size = rj.last() - rj.first() + 1;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int j=r_j.first(); j<=r_j.last(); j++) {
    mons[c++] = new Location(NULL, d, i, j, k);
  }
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n   - rack's name (it is printed before the values),
*  \param d   - Domain on which the rack is created,
*  \param rk  - range of i's (minimum and maximum)
*  \param i,j - logical coordinates defining racks's position.
*
*  \note This constructor creates the rack in "k" direction.                     
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d, 
           const int i, const int j, const Range<int> rk) 
 : name(n), r_i(i,i), r_j(j,j), r_k(rk) {

  int size = rk.last() - rk.first() + 1;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int k=r_k.first(); k<=r_k.last(); k++) {
    mons[c++] = new Location(NULL, d, i, j, k);
  }
}

/******************************************************************************/
void Rack::print(const Scalar & phi) {

  if(name)
    boil::oout << "Rack:" << name << boil::endl;

  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) 
    for(int j=r_j.first(); j<=r_j.last(); j++) 
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        (mons[c++])->print(phi);
      }
}

/******************************************************************************/
void Rack::print(const Vector & u, const Comp & m) {

  if(name)
    boil::oout << "Rack:" << name << boil::endl;

  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) 
    for(int j=r_j.first(); j<=r_j.last(); j++) 
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        (mons[c++])->print(u,m);
      }
}

/*-----------------------------------------------------------------------------+
 '$Id: rack.cpp,v 1.8 2009/05/30 07:45:12 niceno Exp $'/
+-----------------------------------------------------------------------------*/
