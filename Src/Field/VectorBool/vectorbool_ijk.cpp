#include "vectorbool.h"

/******************************************************************************/
int VectorBool::i(const Comp & m, const real x) const {
  if( x<xn(m,si(m)))   return si(m);
  if( x>xn(m,ei(m)+1)) return ei(m);
  for(int i_=0; i_<ni(m); i_++)
    if(x >= xn(m,i_) && x <= xn(m,i_+1)) return i_;
  return -1;
}

/******************************************************************************/
int VectorBool::j(const Comp & m, const real y) const {
  if( y<yn(m,sj(m)))   return sj(m);
  if( y>yn(m,ej(m)+1)) return ej(m);
  for(int j_=0; j_<nj(m); j_++)
    if(y >= yn(m,j_) && y <= yn(m,j_+1)) return j_;
  return -1;
}

/******************************************************************************/
int VectorBool::k(const Comp & m, const real z) const {
  if( z<zn(m,sk(m)))   return sk(m);
  if( z>zn(m,ek(m)+1)) return ek(m);
  for(int k_=0; k_<nk(m); k_++)
    if(z >= zn(m,k_) && z <= zn(m,k_+1)) return k_;
  return -1;
}

/******************************************************************************/
int VectorBool::im(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m);
  if( x>xn(m,ei(m)+1)) return ei(m);
  for(int i_=0; i_<ni(m); i_++)
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int VectorBool::jm(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m);
  if( y>yn(m,ej(m)+1)) return ej(m);
  for(int j_=0; j_<nj(m); j_++)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int VectorBool::km(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m);
  if( z>zn(m,ek(m)+1)) return ek(m);
  for(int k_=0; k_<nk(m); k_++)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int VectorBool::ip(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m);
  if( x>xn(m,ei(m)+1)) return ei(m);
  for(int i_=ni(m)-1; i_>=0; i_--) 
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int VectorBool::jp(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m);
  if( y>yn(m,ej(m)+1)) return ej(m);
  for(int j_=nj(m)-1; j_>=0; j_--)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int VectorBool::kp(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m);
  if( z>zn(m,ek(m)+1)) return ek(m);
  for(int k_=nk(m)-1; k_>=0; k_--)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int VectorBool::aim(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m)-1;
  if( x>xn(m,ei(m)+1)) return ei(m)+1;
  for(int i_=0; i_<ni(m); i_++)
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int VectorBool::ajm(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m)-1;
  if( y>yn(m,ej(m)+1)) return ej(m)+1;
  for(int j_=0; j_<nj(m); j_++)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int VectorBool::akm(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m)-1;
  if( z>zn(m,ek(m)+1)) return ek(m)+1;
  for(int k_=0; k_<nk(m); k_++)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int VectorBool::aip(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m)-1;
  if( x>xn(m,ei(m)+1)) return ei(m)+1;
  for(int i_=ni(m)-1; i_>=0; i_--)
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int VectorBool::ajp(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m)-1;
  if( y>yn(m,ej(m)+1)) return ej(m)+1;
  for(int j_=nj(m)-1; j_>=0; j_--)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int VectorBool::akp(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m)-1;
  if( z>zn(m,ek(m)+1)) return ek(m)+1;
  for(int k_=nk(m)-1; k_>=0; k_--)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/*-----------------------------------------------------------------------------+
 '$Id: vectorbool_ijk.cpp,v 1.1 2014/02/04 08:20:58 sato Exp $'/
+-----------------------------------------------------------------------------*/
