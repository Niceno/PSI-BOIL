#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::bdphiface_check(const Vector & vec
                      , const Scalar & sca) {
/***************************************************************************//**
*  \brief Check boundary condition for face variables.
*******************************************************************************/

  /* i-min or i-max */
  int i1,i2;
  Comp m=Comp::i();
  for_vmjk(vec,m,j,k){
    i1=vec.si(m);
    i2=vec.ei(m);
    if(vec[m][i1][j][k]!=vec[m][i2][j][k])
      std::cout<<" i-direction "<<j<<" "<<k<<" "<<vec[m][i1][j][k]<<" "
               <<vec[m][i2][j][k]<<" "<<vec[m][i1][j][k]-vec[m][i2][j][k]<<"\n";
  }

  /* j-min or j-max */
  int j1,j2;
  m=Comp::j();
  for_vmik(vec,m,i,k){
    j1=vec.sj(m);
    j2=vec.ej(m);
    if(vec[m][i][j1][k]!=vec[m][i][j2][k])
      std::cout<<" j-direction "<<i<<" "<<k<<" "<<vec[m][i][j1][k]<<" "
               <<vec[m][i][j2][k]<<" "<<vec[m][i][j1][k]-vec[m][i][j2][k]<<"\n";
  }

  /* k-min or k-max */
  int k1,k2;
  m=Comp::k();
  for_vmij(vec,m,i,j){
    k1=vec.sk(m);
    k2=vec.ek(m);
    if(vec[m][i][j][k1]!=vec[m][i][j][k2])
      std::cout<<" k-direction "<<i<<" "<<j<<" "<<vec[m][i][j][k1]<<" "
               <<vec[m][i][j][k2]<<" "<<vec[m][i][j][k1]-vec[m][i][j][k2]<<"\n";
  }
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_bdphiface_check.cpp,v 1.1 2010/06/02 09:26:30 sato Exp $'/
+-----------------------------------------------------------------------------*/
