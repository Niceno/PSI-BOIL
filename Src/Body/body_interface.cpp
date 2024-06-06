#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Field/ScalarBool/scalarbool.h"
#include "../Field/VectorBool/vectorbool.h"
#include "../Plot/plot.h"

/******************************************************************************/
real Body::fV(const int i, const int j, const int k) const {

  return (*sca)[i][j][k];
}

/******************************************************************************/
real Body::fV(const Comp & m, const int i, const int j, const int k) const {

  return (*vec)[m][i][j][k];
}

/******************************************************************************/
real Body::fV(const int cc) const {

  assert(cc < cells[3].size());

  return cells[3][cc].fV();
}

/******************************************************************************/
real Body::fV(const Comp & m, const int cc) const {

  assert(cc < cells[~m].size());

  return cells[~m][cc].fV();
}

/******************************************************************************/
bool Body::cut(int i, int j, int k) const {
/*-----------------------------------------+
|  return true if cell is cut, but active  |
+-----------------------------------------*/
  const int cc = index[3][i][j][k];
  if(cc == -1){
    return false;
  } else {
    if(fV(i,j,k) < 0.5){
      return false;
    } else {
      return true;
    }
  }
  //return (fV(i,j,k) > 0.5 && fV(i,j,k) < 1.0);
}

/******************************************************************************/
bool Body::cut_p(int i, int j, int k) const {
/*-----------------------------------------+
|  return true if cell is cut, but active  |
+-----------------------------------------*/
  return (fV(i,j,k) > 0.0 && fV(i,j,k) < 1.0);
}

/******************************************************************************/
bool Body::cut(const Comp & m, int i, int j, int k) const {
/*-----------------------------------------+
|  return true if cell is cut, but active  |
+-----------------------------------------*/
  const int cc = index[~m][i][j][k];
  if(cc == -1){ 
    return false;
  } else {
    if(fV(m,i,j,k) < 0.5){
      return false;
    } else {
      return true;
    }
  }
  //return (fV(m,i,j,k) > 0.5 && fV(m,i,j,k) < 1.0);
}

/******************************************************************************/
bool Body::cut(const int cc) const {
/*-----------------------------------------+
|  return true if cell is cut, but active  |
+-----------------------------------------*/
  assert(cc < cells[3].size());
  if(cc == -1){
    return false;
  } else {
    if(cells[3][cc].fV() < 0.5){
      return false;
    } else {
      return true;
    }
  }
  //return cells[3][cc].fV() > 0.5 && cells[3][cc].fV() < 1.0;
}

/******************************************************************************/
bool Body::cut(const Comp & m, int cc) const {
/*-----------------------------------------+
|  return true if cell is cut, but active  |
+-----------------------------------------*/
  assert(cc < cells[~m].size());
  if(cc == -1){
    return false;
  } else {
    if(cells[~m][cc].fV() < 0.5){
      return false;
    } else {
      return true;
    }
  }
  //return cells[~m][cc].fV() > 0.5 && cells[~m][cc].fV() < 1.0;
}

/******************************************************************************/
bool Body::off(int i, int j, int k) const {
/*--------------------------------------------------------------+
|  return true if cell is inactive, i.e. it is inside the body  |
+--------------------------------------------------------------*/

  return (fV(i,j,k) <= 0.5 + boil::pico);
}

/******************************************************************************/
bool Body::off_p(int i, int j, int k) const {
/*-----------------------------------------------------------------------+
|  return true if pressure cell is inactive, i.e. it is inside the body  |
+-----------------------------------------------------------------------*/

  //return (fV(i,j,k) < boil::pico);
  return (fV(i,j,k) == 0.0);
}

/******************************************************************************/
bool Body::off(const Comp & m, int i, int j, int k) const {
/*--------------------------------------------------------------+
|  return true if cell is inactive, i.e. it is inside the body  |
+--------------------------------------------------------------*/
#if 0
  return (fV(m,i,j,k) <= 0.5 + boil::pico);
#endif
#if 0
  int nn = ~m - 1;
  int ii = int(-0.5*real(nn-abs(nn)));
  int jj = ~m % 2;
  int kk = int( 0.5*real(nn+abs(nn)));
  std::cout<<m<<" "<<i<<" "<<j<<" "<<k<<"\n";
  return (0.5*(fV(i,j,k)+fV(i-ii,j-jj,k-kk)) <= 0.5 + boil::pico);
#endif
#if 1
  return (*vecoff)[m][i][j][k];
#endif
}

/******************************************************************************/
bool Body::on(int i, int j, int k) const {
/*----------------------------------------------------------+
|  return true if cell is active, i.e. it is ins the fluid  |
+----------------------------------------------------------*/

  return !off(i,j,k);
}

/******************************************************************************/
bool Body::on_p(int i, int j, int k) const {
/*-------------------------------------------------------------------+
|  return true if pressure cell is active, i.e. it is ins the fluid  |
+-------------------------------------------------------------------*/

  return !off_p(i,j,k);
}

/******************************************************************************/
bool Body::on(const Comp & m, int i, int j, int k) const {
/*----------------------------------------------------------+
|  return true if cell is active, i.e. it is ins the fluid  |
+----------------------------------------------------------*/

  return !off(m,i,j,k);
}

/******************************************************************************/
real Body::dSx(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1) return 0.0; // is this the best practice?
  else         return polys[cc].area_x();
}

/******************************************************************************/
real Body::dSy(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1) return 0.0; // is this the best practice?
  else         return polys[cc].area_y();
}

/******************************************************************************/
real Body::dSz(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1) return 0.0; // is this the best practice?
  else         return polys[cc].area_z();
}

/******************************************************************************/
real Body::fSw(const int cc) const {return cells[3][cc].fS(0,0);}
real Body::fSe(const int cc) const {return cells[3][cc].fS(0,1);}
real Body::fSs(const int cc) const {return cells[3][cc].fS(1,0);}
real Body::fSn(const int cc) const {return cells[3][cc].fS(1,1);}
real Body::fSb(const int cc) const {return cells[3][cc].fS(2,0);}
real Body::fSt(const int cc) const {return cells[3][cc].fS(2,1);}

/******************************************************************************/
real Body::fPmmm(const int cc) const {return cells[3][cc].fP(0,0,0);}
real Body::fPpmm(const int cc) const {return cells[3][cc].fP(1,0,0);}
real Body::fPmpm(const int cc) const {return cells[3][cc].fP(0,1,0);}
real Body::fPppm(const int cc) const {return cells[3][cc].fP(1,1,0);}
real Body::fPmmp(const int cc) const {return cells[3][cc].fP(0,0,1);}
real Body::fPpmp(const int cc) const {return cells[3][cc].fP(1,0,1);}
real Body::fPmpp(const int cc) const {return cells[3][cc].fP(0,1,1);}
real Body::fPppp(const int cc) const {return cells[3][cc].fP(1,1,1);}

/******************************************************************************/
real Body::fE000(const int cc) const {return cells[3][cc].fE(0,0,0);}
real Body::fE010(const int cc) const {return cells[3][cc].fE(0,1,0);}
real Body::fE001(const int cc) const {return cells[3][cc].fE(0,0,1);}
real Body::fE011(const int cc) const {return cells[3][cc].fE(0,1,1);}
real Body::fE100(const int cc) const {return cells[3][cc].fE(1,0,0);}
real Body::fE110(const int cc) const {return cells[3][cc].fE(1,1,0);}
real Body::fE101(const int cc) const {return cells[3][cc].fE(1,0,1);}
real Body::fE111(const int cc) const {return cells[3][cc].fE(1,1,1);}
real Body::fE200(const int cc) const {return cells[3][cc].fE(2,0,0);}
real Body::fE210(const int cc) const {return cells[3][cc].fE(2,1,0);}
real Body::fE201(const int cc) const {return cells[3][cc].fE(2,0,1);}
real Body::fE211(const int cc) const {return cells[3][cc].fE(2,1,1);}

/******************************************************************************/
int Body::iacpt(const int cc) const {return cells[3][cc].iacpt();}
int Body::jacpt(const int cc) const {return cells[3][cc].jacpt();}
int Body::kacpt(const int cc) const {return cells[3][cc].kacpt();}

/******************************************************************************/
real Body::fSw(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fS(0,0);
}

/******************************************************************************/
real Body::fSe(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fS(0,1);
}

/******************************************************************************/
real Body::fSs(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fS(1,0);
}

/******************************************************************************/
real Body::fSn(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fS(1,1);
}

/******************************************************************************/
real Body::fSb(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fS(2,0);
}

/******************************************************************************/
real Body::fSt(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fS(2,1);
}

/******************************************************************************/
real Body::fPmmm(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==1.0){return 1.0;}
    else {return 0.0;}
  else
    return cells[3][cc].fP(0,0,0);
}

/******************************************************************************/
real Body::fPpmm(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(1,0,0);
}

/******************************************************************************/
real Body::fPmpm(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(0,1,0);
}

/******************************************************************************/
real Body::fPppm(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(1,1,0);
}

/******************************************************************************/
real Body::fPmmp(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(0,0,1);
}

/******************************************************************************/
real Body::fPpmp(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(1,0,1);
}

/******************************************************************************/
real Body::fPmpp(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1) 
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(0,1,1);
}   

/******************************************************************************/
real Body::fPppp(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fP(1,1,1);
}

/******************************************************************************/
real Body::fE000(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(0,0,0);
}

/******************************************************************************/
real Body::fE010(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(0,1,0);
}

/******************************************************************************/
real Body::fE001(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(0,0,1);
}

/******************************************************************************/
real Body::fE011(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(0,1,1);
}

/******************************************************************************/
real Body::fE100(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(1,0,0);
}

/******************************************************************************/
real Body::fE110(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(1,1,0);
}

/******************************************************************************/
real Body::fE101(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(1,0,1);
}

/******************************************************************************/
real Body::fE111(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(1,1,1);
}

/******************************************************************************/
real Body::fE200(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(2,0,0);
}

/******************************************************************************/
real Body::fE210(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(2,1,0);
}

/******************************************************************************/
real Body::fE201(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(2,0,1);
}

/******************************************************************************/
real Body::fE211(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[3][cc].fE(2,1,1);
}

/******************************************************************************/
int Body::iacpt(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return -1;
  else
    return cells[3][cc].iacpt();
}

/******************************************************************************/
int Body::jacpt(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return -1;
  else
    return cells[3][cc].jacpt();
}
/******************************************************************************/
int Body::kacpt(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return -1;
  else
    return cells[3][cc].kacpt();
}

/******************************************************************************/
real Body::fdxw(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return 1.0;
  else
    return cells[3][cc].fdxw();
}

/******************************************************************************/
real Body::fdxe(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return 1.0;
  else
    return cells[3][cc].fdxe();
}

/******************************************************************************/
real Body::fdys(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return 1.0;
  else
    return cells[3][cc].fdys();
}

/******************************************************************************/
real Body::fdyn(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return 1.0;
  else
    return cells[3][cc].fdyn();
}

/******************************************************************************/
real Body::fdzb(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return 1.0;
  else
    return cells[3][cc].fdzb();
}

/******************************************************************************/
real Body::fdzt(const int i, const int j, const int k) const {
  const int cc = index[3][i][j][k];
  if(cc == -1)
    return 1.0;
  else
    return cells[3][cc].fdzt();
}

/******************************************************************************/
real Body::fSw(const Comp & m, const int cc) const 
                                                 {return cells[~m][cc].fS(0,0);}
real Body::fSe(const Comp & m, const int cc) const 
                                                 {return cells[~m][cc].fS(0,1);}
real Body::fSs(const Comp & m, const int cc) const 
                                                 {return cells[~m][cc].fS(1,0);}
real Body::fSn(const Comp & m, const int cc) const 
                                                 {return cells[~m][cc].fS(1,1);}
real Body::fSb(const Comp & m, const int cc) const 
                                                 {return cells[~m][cc].fS(2,0);}
real Body::fSt(const Comp & m, const int cc) const 
                                                 {return cells[~m][cc].fS(2,1);}

/******************************************************************************/
real Body::fSw(const Comp & m, const int i, const int j, const int k) const {
  const int cc = index[~m][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(m,i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[~m][cc].fS(0,0);
}

/******************************************************************************/
real Body::fSe(const Comp & m, const int i, const int j, const int k) const {
  const int cc = index[~m][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(m,i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[~m][cc].fS(0,1);
}

/******************************************************************************/
real Body::fSs(const Comp & m, const int i, const int j, const int k) const {
  const int cc = index[~m][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(m,i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[~m][cc].fS(1,0);
}

/******************************************************************************/
real Body::fSn(const Comp & m, const int i, const int j, const int k) const {
  const int cc = index[~m][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(m,i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[~m][cc].fS(1,1);
}

/******************************************************************************/
real Body::fSb(const Comp & m, const int i, const int j, const int k) const {
  const int cc = index[~m][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(m,i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[~m][cc].fS(2,0);
}

/******************************************************************************/
real Body::fSt(const Comp & m, const int i, const int j, const int k) const {
  const int cc = index[~m][i][j][k];
  if(cc == -1)
    //return 1.0;
    if(fV(m,i,j,k)==0.0){return 0.0;}
    else {return 1.0;}
  else
    return cells[~m][cc].fS(2,1);
}

/******************************************************************************/
real Body::fdxw(const Comp & m, const int cc) const 
                                                  {return cells[~m][cc].fdxw();}
real Body::fdxe(const Comp & m, const int cc) const 
                                                  {return cells[~m][cc].fdxe();}
real Body::fdys(const Comp & m, const int cc) const 
                                                  {return cells[~m][cc].fdys();}
real Body::fdyn(const Comp & m, const int cc) const 
                                                  {return cells[~m][cc].fdyn();}
real Body::fdzb(const Comp & m, const int cc) const 
                                                  {return cells[~m][cc].fdzb();}
real Body::fdzt(const Comp & m, const int cc) const 
                                                  {return cells[~m][cc].fdzt();}

/******************************************************************************/
real Body::fdxw(const int cc) const {return cells[3][cc].fdxw();}
real Body::fdxe(const int cc) const {return cells[3][cc].fdxe();}
real Body::fdys(const int cc) const {return cells[3][cc].fdys();}
real Body::fdyn(const int cc) const {return cells[3][cc].fdyn();}
real Body::fdzb(const int cc) const {return cells[3][cc].fdzb();}
real Body::fdzt(const int cc) const {return cells[3][cc].fdzt();}

/******************************************************************************/
int Body::nccells() const {

  //if(cells.size() > 0) return cells[3].size();
  if(cells.size() > 0) {return nccells_in[3];}
  else                 {return 0;}
}

/******************************************************************************/
int Body::ncall() const {
  if(cells.size() > 0) {return ncall_;}
  else                 {return 0;}
}

/******************************************************************************/
int Body::nccells(const Comp & m) const {

  //if(cells.size() > 0) return cells[~m].size();
  if(cells.size() > 0) {return nccells_in[~m];}
  else                 {return 0;}
}

/******************************************************************************/
void Body::ijk(int cc, int * i, int * j, int * k) const  {
 
  assert(cc < cells[3].size());

  return cells[3][cc].ijk(i,j,k);
}

/******************************************************************************/
void Body::ijk(Comp & m, int cc, int * i, int * j, int * k) const  {
 
  assert(cc < cells[~m].size());

  return cells[~m][cc].ijk(i,j,k);
}

/******************************************************************************/
real Body::dist(const int i, const int j, const int k) const {

  return (*bdist)[i][j][k];
}

/******************************************************************************/
real Body::nwx(const int i, const int j, const int k) const {

  return (*ux)[i][j][k];
}

/******************************************************************************/
real Body::nwy(const int i, const int j, const int k) const {

  return (*uy)[i][j][k];
}

/******************************************************************************/
real Body::nwz(const int i, const int j, const int k) const {

  return (*uz)[i][j][k];
}
