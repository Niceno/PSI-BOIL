#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::calculate_node_temperature(const Scalar * diff_eddy) {
/***************************************************************************//*** 
*  \brief calculate node temperature in the domain
*  bndtpr = boil::unreal indicates averaging across the interface (= banned)
*******************************************************************************/
  Comp m;
 
  m = Comp::i();
  for(int i=si(); i<=ei()+1; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()  ; k++) {
    real temp = boil::unreal;
    if(!Interface(+1,m,i-1,j,k)) {
#if 0
      temp = tpr[i][j][k]*tpr.dxc(i-1) + tpr[i-1][j][k]*tpr.dxc(i);
      temp /= (tpr.dxc(i) + tpr.dxc(i-1));
#else
      real len_1 = phi.dxw(i) - 0.5*phi.dxc(i);
      real len_2 = 0.5*phi.dxc(i);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i-1][j][k];

      real lam_1,lam_2;

      if(dom->ibody().off(i,j,k)) {
        lam_1 = solid()->lambda(i,j,k);
      } else {
        if(clr[i][j][k]<clrsurf) {
          lam_1=lambdav;
          if(diff_eddy) lam_1 += (*diff_eddy)[i][j][k]*cpv/rhov/turbP; 
        } else {
          lam_1=lambdal;
          if(diff_eddy) lam_1 += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
        }
      }
  
      if(dom->ibody().off(i-1,j,k)) {
        lam_2 = solid()->lambda(i-1,j,k);
      } else {
        if(clr[i-1][j][k]<clrsurf) {
          lam_2=lambdav;
          if(diff_eddy) lam_2 += (*diff_eddy)[i-1][j][k]*cpv/rhov/turbP;
        } else {
          lam_2=lambdal;
          if(diff_eddy) lam_2 += (*diff_eddy)[i-1][j][k]*cpl/rhol/turbP;
        }
      }

      temp = temperature_node(len_1, lam_1, tpr_1, len_2, lam_2, tpr_2);
#endif
    }
    bndtpr[m][i][j][k] = temp;
  }

  m = Comp::j();
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()+1; j++)
  for(int k=sk(); k<=ek()  ; k++) {
    real temp = boil::unreal;
    if(!Interface(+1,m,i,j-1,k)) {
#if 0
      temp = tpr[i][j][k]*tpr.dyc(j-1) + tpr[i][j-1][k]*tpr.dyc(j);
      temp /= (tpr.dyc(j) + tpr.dyc(j-1));
#else
      real len_1 = phi.dys(j) - 0.5*phi.dyc(j);
      real len_2 = 0.5*phi.dyc(j);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j-1][k];

      real lam_1,lam_2;

      if(dom->ibody().off(i,j,k)) {
        lam_1 = solid()->lambda(i,j,k);
      } else {
        if(clr[i][j][k]<clrsurf) {
          lam_1=lambdav;
          if(diff_eddy) lam_1 += (*diff_eddy)[i][j][k]*cpv/rhov/turbP; 
        } else {
          lam_1=lambdal;
          if(diff_eddy) lam_1 += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
        }
      }
  
      if(dom->ibody().off(i,j-1,k)) {
        lam_2 = solid()->lambda(i,j-1,k);
      } else {
        if(clr[i][j-1][k]<clrsurf) {
          lam_2=lambdav;
          if(diff_eddy) lam_2 += (*diff_eddy)[i][j-1][k]*cpv/rhov/turbP;
        } else {
          lam_2=lambdal;
          if(diff_eddy) lam_2 += (*diff_eddy)[i][j-1][k]*cpl/rhol/turbP;
        }
      }

      temp = temperature_node(len_1, lam_1, tpr_1, len_2, lam_2, tpr_2);
#endif
    }
    bndtpr[m][i][j][k] = temp;
  }

  m = Comp::k();
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()+1; k++) {
    real temp = boil::unreal;
    if(!Interface(+1,m,i,j,k-1)) {
#if 0
      temp = tpr[i][j][k]*tpr.dzc(k-1) + tpr[i][j][k-1]*tpr.dzc(k);
      temp /= (tpr.dzc(k) + tpr.dzc(k-1));
#else
      real len_1 = phi.dzb(k) - 0.5*phi.dzc(k);
      real len_2 = 0.5*phi.dzc(k);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j][k-1];

      real lam_1,lam_2;

      if(dom->ibody().off(i,j,k)) {
        lam_1 = solid()->lambda(i,j,k);
      } else {
        if(clr[i][j][k]<clrsurf) {
          lam_1=lambdav;
          if(diff_eddy) lam_1 += (*diff_eddy)[i][j][k]*cpv/rhov/turbP; 
        } else {
          lam_1=lambdal;
          if(diff_eddy) lam_1 += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
        }
      }
  
      if(dom->ibody().off(i,j,k-1)) {
        lam_2 = solid()->lambda(i,j,k-1);
      } else {
        if(clr[i][j][k-1]<clrsurf) {
          lam_2=lambdav;
          if(diff_eddy) lam_2 += (*diff_eddy)[i][j][k-1]*cpv/rhov/turbP;
        } else {
          lam_2=lambdal;
          if(diff_eddy) lam_2 += (*diff_eddy)[i][j][k-1]*cpl/rhol/turbP;
        }
      }

      temp = temperature_node(len_1, lam_1, tpr_1, len_2, lam_2, tpr_2);
#endif
    }
    bndtpr[m][i][j][k] = temp;
  }

  bndtpr.exchange_all();
}
