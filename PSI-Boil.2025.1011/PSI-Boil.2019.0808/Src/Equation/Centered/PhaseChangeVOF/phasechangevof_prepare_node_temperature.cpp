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
      real len_1 = 0.5*phi.dxc(i);
      real len_2 = 0.5*phi.dxc(i-1);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i-1][j][k];

      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i-1,j,k,diff_eddy);

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
      real len_1 = 0.5*phi.dyc(j);
      real len_2 = 0.5*phi.dyc(j-1);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j-1][k];

      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j-1,k,diff_eddy);

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
      real len_1 = 0.5*phi.dzc(k);
      real len_2 = 0.5*phi.dzc(k-1);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j][k-1];

      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j,k-1,diff_eddy);

      temp = temperature_node(len_1, lam_1, tpr_1, len_2, lam_2, tpr_2);
#endif
    }
    bndtpr[m][i][j][k] = temp;
  }

#if 0
  for_vijk(clr,i,j,k) {
    //if(adens[i][j][k]>boil::atto) {
    if(i<10&&j<10&&k<10 && i>boil::BW&&j>boil::BW&&k>boil::BW) {
      boil::oout<<i<<" "<<j<<" "<<k<<" "<<adens[i][j][k]<<" | "<<tpr[i][j][k-1]<<" "<<tpr[i][j][k]<<" "<<tpr[i][j][k+1]<<" "
                <<"| "<<bndtpr[Comp::w()][i][j][k]<<" "<<bndtpr[Comp::w()][i][j][k+1]<<" | "
                <<(lambda(i,j,k-1,diff_eddy)*tpr[i][j][k-1]*0.5*phi.dzc(k)+lambda(i,j,k,diff_eddy)*tpr[i][j][k]*0.5*phi.dzc(k-1))/(lambda(i,j,k-1,diff_eddy)*0.5*phi.dzc(k)+lambda(i,j,k,diff_eddy)*0.5*phi.dzc(k-1))<<" "
                <<(lambda(i,j,k+1,diff_eddy)*tpr[i][j][k+1]*0.5*phi.dzc(k)+lambda(i,j,k,diff_eddy)*tpr[i][j][k]*0.5*phi.dzc(k+1))/(lambda(i,j,k+1,diff_eddy)*0.5*phi.dzc(k)+lambda(i,j,k,diff_eddy)*0.5*phi.dzc(k+1))<<" | "
                <<lambda(i,j,k+1,diff_eddy)*(tpr[i][j][k+1]-bndtpr[Comp::w()][i][j][k+1])/(0.5*phi.dzc(k+1))
                 -lambda(i,j,k,diff_eddy)*(bndtpr[Comp::w()][i][j][k+1]-tpr[i][j][k])/(0.5*phi.dzc(k))
                <<" "
                <<lambda(i,j,k,diff_eddy)*(tpr[i][j][k]-bndtpr[Comp::w()][i][j][k])/(0.5*phi.dzc(k))
                 -lambda(i,j,k-1,diff_eddy)*(bndtpr[Comp::w()][i][j][k]-tpr[i][j][k-1])/(0.5*phi.dzc(k-1))
                <<" | "<<lambda(i,j,k+1,diff_eddy)*(tpr[i][j][k+1]-bndtpr[Comp::w()][i][j][k+1])/(0.5*phi.dzc(k+1))<<" "<<lambda(i,j,k,diff_eddy)*(bndtpr[Comp::w()][i][j][k+1]-tpr[i][j][k])/(0.5*phi.dzc(k))
<<" "<<lambda(i,j,k,diff_eddy)*(tpr[i][j][k]-bndtpr[Comp::w()][i][j][k])/(0.5*phi.dzc(k))<<" "<<lambda(i,j,k-1,diff_eddy)*(bndtpr[Comp::w()][i][j][k]-tpr[i][j][k-1])/(0.5*phi.dzc(k-1))
                <<" | "<<lambda(i,j,k-1,diff_eddy)<<" "<<lambda(i,j,k,diff_eddy)<<" "<<lambda(i,j,k+1,diff_eddy)<<boil::endl;
    }
  }
  exit(0);
#endif

  bndtpr.exchange_all();
}
