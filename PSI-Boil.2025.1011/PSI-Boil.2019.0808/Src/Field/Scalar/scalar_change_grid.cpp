#include "scalar.h"

/******************************************************************************/
void Scalar::change_grid(const Scalar & s_org) {
/*---------------------+
|  get val from s_org  |
+---------------------*/
  for_ijk(i,j,k){
    real x_new = xc(i);
    real y_new = yc(j);
    real z_new = zc(k);
    int i_org=-1, j_org=-1, k_org=-1;
    for_avi(s_org,ii){if((s_org.xn(ii)<=x_new)&&(x_new<s_org.xn(ii+1))) {i_org = ii; break;}}
    for_avj(s_org,jj){if((s_org.yn(jj)<=y_new)&&(y_new<s_org.yn(jj+1))) {j_org = jj; break;}}
    for_avk(s_org,kk){if((s_org.zn(kk)<=z_new)&&(z_new<s_org.zn(kk+1))) {k_org = kk; break;}}
    assert(i_org>=0);
    assert(j_org>=0);
    assert(k_org>=0);
    val[i][j][k]=s_org[i_org][j_org][k_org];
  }
}
