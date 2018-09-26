#include "scalarbool.h"

/******************************************************************************/
int ScalarBool::i(const real x) const {
  if( x<xn(si()))   return si();
  if( x>xn(ei()+1)) return ei();
  for(int i_=0; i_<ni(); i_++) if(x >= xn(i_) && x <= xn(i_+1)) return i_;
  return -1;
}

/******************************************************************************/
int ScalarBool::j(const real y) const {
  if( y<yn(sj()))   return sj();
  if( y>yn(ej()+1)) return ej();
  for(int j_=0; j_<nj(); j_++) if(y >= yn(j_) && y <= yn(j_+1)) return j_;
  return -1;
}

/******************************************************************************/
int ScalarBool::k(const real z) const {
  if( z<zn(sk()))   return sk();
  if( z>zn(ek()+1)) return ek();
  for(int k_=0; k_<nk(); k_++) if(z >= zn(k_) && z <= zn(k_+1)) return k_;
  return -1;
}

/******************************************************************************/
int ScalarBool::im(const real x, const real t) const{
  if( x+t<xn(si()))   return si();
  if( x-t>xn(ei()+1)) return ei();
  for(int i_=1; i_<ni(); i_++) if(x >= xn(i_)-t && x <= xn(i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int ScalarBool::jm(const real y, const real t) const{
  if( y+t<yn(sj()))   return sj();
  if( y-t>yn(ej()+1)) return ej();
  for(int j_=1; j_<nj(); j_++) if(y >= yn(j_)-t && y <= yn(j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int ScalarBool::km(const real z, const real t) const{
  if( z+t<zn(sk()))   return sk();
  if( z-t>zn(ek()+1)) return ek();
  for(int k_=1; k_<nk(); k_++) if(z >= zn(k_)-t && z <= zn(k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int ScalarBool::ip(const real x, const real t) const {
  if( x+t<xn(si()))   return si();
  if( x-t>xn(ei()+1)) return ei();
  for(int i_=ni()-1; i_>=0; i_--)if(x >= xn(i_)-t && x <= xn(i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int ScalarBool::jp(const real y, const real t) const {
  if( y+t<yn(sj()))   return sj();
  if( y-t>yn(ej()+1)) return ej();
  for(int j_=nj()-1; j_>=0; j_--)if(y >= yn(j_)-t && y <= yn(j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int ScalarBool::kp(const real z, const real t) const {
  if( z+t<zn(sk()))   return sk();
  if( z-t>zn(ek()+1)) return ek();
  for(int k_=nk()-1; k_>=0; k_--)if(z >= zn(k_)-t && z <= zn(k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int ScalarBool::aim(const real x, const real t) const {
  if( x+t<xn(si()-1)) return si()-1;
  if( x-t>xn(ei()+2)) return ei()+1;
  for(int i_=0; i_<ni(); i_++)if(x >= xn(i_)-t && x <= xn(i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int ScalarBool::ajm(const real y, const real t) const {
  if( y+t<yn(sj()-1)) return sj()-1;
  if( y-t>yn(ej()+2)) return ej()+1;
  for(int j_=0; j_<nj(); j_++)if(y >= yn(j_)-t && y <= yn(j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int ScalarBool::akm(const real z, const real t) const {
  if( z+t<zn(sk()-1)) return sk()-1;
  if( z-t>zn(ek()+2)) return ek()+1;
  for(int k_=0; k_<nk(); k_++)if(z >= zn(k_)-t && z <= zn(k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int ScalarBool::aip(const real x, const real t) const {
  if( x+t<xn(si()-1)) return si()-1;
  if( x-t>xn(ei()+2)) return ei()+1;
  for(int i_=ni()-1; i_>=0; i_--)if(x >= xn(i_)-t && x <= xn(i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int ScalarBool::ajp(const real y, const real t) const {
  if( y+t<yn(sj()-1)) return sj()-1;
  if( y-t>yn(ej()+2)) return ej()+1;
  for(int j_=nj()-1; j_>=0; j_--)if(y >= yn(j_)-t && y <= yn(j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int ScalarBool::akp(const real z, const real t) const {
  if( z+t<zn(sk()-1)) return sk()-1;
  if( z-t>zn(ek()+2)) return ek()+1;
  for(int k_=nk()-1; k_>=0; k_--)if(z >= zn(k_)-t && z <= zn(k_+1)+t) return k_;
  return -1;
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarbool_ijk.cpp,v 1.2 2017/08/09 08:33:15 sato Exp $'/
+-----------------------------------------------------------------------------*/
