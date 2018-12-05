#include "finescalar.h"

/******************************************************************************/
void FineScalar::eval_face() {
/***************************************************************************//**
*  \brief Evaluate the marker function at faces.
*******************************************************************************/

  /* simplified way, assuming fs was already evaluated */
  Comp m;

  /*------+
  | x-dir |
  +------*/
  m = Comp::i();

  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()  ; j++)
  for(int k=phi->sk(); k<=phi->ek()  ; k++) {

    real fsval = (*fs)[m][i][j][k];
    real phim = (*phi)[i-1][j][k];
    real phip = (*phi)[i  ][j][k];
    
    if(!boil::realistic(fsval)) {
      value(i,j,k,w()) = phim>phisurf;
    } else {
      real edgepos = phi->xn(i);
      if(fsval>edgepos) {
        value(i,j,k,w()) = phim>phisurf; 
      } else if(fsval<edgepos) {
        value(i,j,k,w()) = phip>phisurf; 
      } else {
        value(i,j,k,w()) = 0.5; 
      }
    }
  }

  /*------+
  | y-dir |
  +------*/
  m = Comp::j();

  for(int i=phi->si(); i<=phi->ei()  ; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()  ; k++) {

    real fsval = (*fs)[m][i][j][k];
    real phim = (*phi)[i][j-1][k];
    real phip = (*phi)[i][j  ][k];
    
    if(!boil::realistic(fsval)) {
      value(i,j,k,s()) = phim>phisurf;
    } else {
      real edgepos = phi->yn(j);
      if(fsval>edgepos) {
        value(i,j,k,s()) = phim>phisurf; 
      } else if(fsval<edgepos) {
        value(i,j,k,s()) = phip>phisurf; 
      } else {
        value(i,j,k,s()) = 0.5; 
      }
    }
  }

  /*------+
  | z-dir |
  +------*/
  m = Comp::k();

  for(int i=phi->si(); i<=phi->ei()  ; i++)
  for(int j=phi->sj(); j<=phi->ej()  ; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    real fsval = (*fs)[m][i][j][k];
    real phim = (*phi)[i][j][k-1];
    real phip = (*phi)[i][j][k  ];
    
    if(!boil::realistic(fsval)) {
      value(i,j,k,b()) = phim>phisurf;
    } else {
      real edgepos = phi->zn(k);
      if(fsval>edgepos) {
        value(i,j,k,b()) = phim>phisurf; 
      } else if(fsval<edgepos) {
        value(i,j,k,b()) = phip>phisurf; 
      } else {
        value(i,j,k,b()) = 0.5; 
      }
    }
  }

  /* correct at boundaries */
  bdcond_face();

}
