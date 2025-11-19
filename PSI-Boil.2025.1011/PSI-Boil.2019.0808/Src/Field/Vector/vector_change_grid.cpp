#include "vector.h"

/******************************************************************************/
void Vector::change_grid(const Vector & v_org) {

  for_m(m){
    for_mijk(m,i,j,k){
      real x_new = xc(m,i);
      real y_new = yc(m,j);
      real z_new = zc(m,k);
      int i_org=-1, j_org=-1, k_org=-1;
      for_avmi(v_org,m,ii){if((v_org.xn(m,ii)<=x_new)&&(x_new<=v_org.xn(m,ii+1))) {i_org = ii; break;}}
      for_avmj(v_org,m,jj){if((v_org.yn(m,jj)<=y_new)&&(y_new<=v_org.yn(m,jj+1))) {j_org = jj; break;}}
      for_avmk(v_org,m,kk){if((v_org.zn(m,kk)<=z_new)&&(z_new<=v_org.zn(m,kk+1))) {k_org = kk; break;}}
      assert(i_org>=0);
      assert(j_org>=0);
      assert(k_org>=0);
      //if (k_org<0){
      //  std::cout<<"m,i,j,k= "<<m<<" "<<i<<" "<<j<<" "<<k<<" "<<z_new<<"\n";
      //  exit(0);
      //}
      vec[m][i][j][k]=v_org[m][i_org][j_org][k_org];
    }
  }
}
