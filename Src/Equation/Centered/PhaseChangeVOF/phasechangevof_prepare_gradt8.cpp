#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::prepare_gradt8() {
/***************************************************************************//*** 
*  \brief prepare gradt in the domain
*  bndtpr = -1000.0 indicates averaging across the interface (= banned)
*******************************************************************************/
  Comp m;
 
  m = Comp::i();
  for(int i=si(); i<=ei()+1; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()  ; k++) {
    real temp = tempnull;
    if(!Interface(i-1,j,k,m,+1)) {
      temp = tpr[i][j][k]*tpr.dxc(i-1) + tpr[i-1][j][k]*tpr.dxc(i);
      temp /= (tpr.dxc(i) + tpr.dxc(i-1));
    }
    bndtpr[m][i][j][k] = temp;
  }

  m = Comp::j();
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()+1; j++)
  for(int k=sk(); k<=ek()  ; k++) {
    real temp = tempnull;
    if(!Interface(i-1,j,k,m,+1)) {
      temp = tpr[i][j][k]*tpr.dyc(j-1) + tpr[i][j-1][k]*tpr.dyc(j);
      temp /= (tpr.dyc(j) + tpr.dyc(j-1));
    }
    bndtpr[m][i][j][k] = temp;
  }

  m = Comp::k();
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()+1; k++) {
    real temp = tempnull;
    if(!Interface(i-1,j,k,m,+1)) {
      temp = tpr[i][j][k]*tpr.dzc(k-1) + tpr[i][j][k-1]*tpr.dzc(k);
      temp /= (tpr.dzc(k) + tpr.dzc(k-1));
    }
    bndtpr[m][i][j][k] = temp;
  }

  bndtpr.exchange_all();
}
