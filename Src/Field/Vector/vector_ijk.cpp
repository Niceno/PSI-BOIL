#include "vector.h"

/******************************************************************************/
int Vector::i(const Comp & m, const real x) const {
  if( x<xn(m,si(m)))   return si(m);
  if( x>xn(m,ei(m)+1)) return ei(m);
  for(int i_=si(m); i_<=ei(m); i_++)
    if(x >= xn(m,i_) && x <= xn(m,i_+1)) return i_;
  return -1;
}

/******************************************************************************/
int Vector::I(const Comp & m, const real x) const {
  int sI = 1;
  int eI = dom->gi()-2;
  if (m==Comp::i()) eI=dom->gi()-1;
  if( x<xn_global(m,sI))   return sI;
  if( x>xn_global(m,eI+1)) return eI;
  for(int I_=sI; I_<=eI; I_++)
    if(x >= xn_global(m,I_) && x <= xn_global(m,I_+1)) return I_;
  return -1;
}

/******************************************************************************/
int Vector::j(const Comp & m, const real y) const {
  if( y<yn(m,sj(m)))   return sj(m);
  if( y>yn(m,ej(m)+1)) return ej(m);
  for(int j_=sj(m); j_<=ej(m); j_++)
    if(y >= yn(m,j_) && y <= yn(m,j_+1)) return j_;
  return -1;
}

/******************************************************************************/
int Vector::J(const Comp & m, const real y) const {
  int sJ = 1;
  int eJ = dom->gj()-2;
  if (m==Comp::j()) eJ=dom->gj()-1;
  if( y<yn_global(m,sJ))   return sJ;
  if( y>yn_global(m,eJ+1)) return eJ;
  for(int J_=sJ; J_<=eJ; J_++)
    if(y >= yn_global(m,J_) && y <= yn_global(m,J_+1)) return J_;
  return -1;
}

/******************************************************************************/
int Vector::k(const Comp & m, const real z) const {
  if( z<zn(m,sk(m)))   return sk(m);
  if( z>zn(m,ek(m)+1)) return ek(m);
  for(int k_=sk(m); k_<=ek(m); k_++)
    if(z >= zn(m,k_) && z <= zn(m,k_+1)) return k_;
  return -1;
}

/******************************************************************************/
int Vector::K(const Comp & m, const real z) const {
  int sK = 1;
  int eK = dom->gk()-2;
  if (m==Comp::k()) eK=dom->gk()-1;
  if( z<zn_global(m,sK))   return sK;
  if( z>zn_global(m,eK+1)) return eK;
  for(int K_=sK; K_<=eK; K_++)
    if(z >= zn_global(m,K_) && z <= zn_global(m,K_+1)) return K_;
  return -1;
}

/******************************************************************************/
int Vector::im(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m);
  if( x>xn(m,ei(m)+1)) return ei(m);
  for(int i_=si(m); i_<=ei(m); i_++)
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int Vector::jm(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m);
  if( y>yn(m,ej(m)+1)) return ej(m);
  for(int j_=sj(m); j_<=ej(m); j_++)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int Vector::km(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m);
  if( z>zn(m,ek(m)+1)) return ek(m);
  for(int k_=sk(m); k_<=ek(m); k_++)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int Vector::ip(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m);
  if( x>xn(m,ei(m)+1)) return ei(m);
  for(int i_=ei(m); i_>=si(m); i_--) 
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int Vector::jp(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m);
  if( y>yn(m,ej(m)+1)) return ej(m);
  for(int j_=ej(m); j_>=sj(m); j_--)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int Vector::kp(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m);
  if( z>zn(m,ek(m)+1)) return ek(m);
  for(int k_=ek(m); k_>=sk(m); k_--)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int Vector::aim(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m)-1;
  if( x>xn(m,ei(m)+1)) return ei(m)+1;
  for(int i_=si(m); i_<=ei(m); i_++)
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int Vector::ajm(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m)-1;
  if( y>yn(m,ej(m)+1)) return ej(m)+1;
  for(int j_=sj(m); j_<=ej(m); j_++)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int Vector::akm(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m)-1;
  if( z>zn(m,ek(m)+1)) return ek(m)+1;
  for(int k_=sk(m); k_<=ek(m); k_++)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/******************************************************************************/
int Vector::aip(const Comp & m, const real x, const real t) const {
  if( x<xn(m,si(m)))   return si(m)-1;
  if( x>xn(m,ei(m)+1)) return ei(m)+1;
  for(int i_=ei(m); i_>=si(m); i_--)
    if(x >= xn(m,i_)-t && x <= xn(m,i_+1)+t) return i_;
  return -1;
}

/******************************************************************************/
int Vector::ajp(const Comp & m, const real y, const real t) const {
  if( y<yn(m,sj(m)))   return sj(m)-1;
  if( y>yn(m,ej(m)+1)) return ej(m)+1;
  for(int j_=ej(m); j_>=sj(m); j_--)
    if(y >= yn(m,j_)-t && y <= yn(m,j_+1)+t) return j_;
  return -1;
}

/******************************************************************************/
int Vector::akp(const Comp & m, const real z, const real t) const {
  if( z<zn(m,sk(m)))   return sk(m)-1;
  if( z>zn(m,ek(m)+1)) return ek(m)+1;
  for(int k_=ek(m); k_>=sk(m); k_--)
    if(z >= zn(m,k_)-t && z <= zn(m,k_+1)+t) return k_;
  return -1;
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_ijk.cpp,v 1.7 2016/03/15 15:37:20 sato Exp $'/
+-----------------------------------------------------------------------------*/
