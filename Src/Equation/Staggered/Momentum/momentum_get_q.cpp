#include "momentum.h" /*spremenil draksler- 29.8. at 10.30*/

/******************************************************************************/
void Momentum::get_q(Scalar * q) {
/*-------------------------------------------------------------------+
|  strain and vorticity tensors are probably not computet properly.  |
|  i think i miss 1/2 in front of each term. check it!               |
+-------------------------------------------------------------------*/

  u.exchange();
  
  for_vijk((*q),i,j,k) {

    real du_dx = (u[Comp::u()][i+1][j][k] - u[Comp::u()][i][j][k]) / u.dxc(i);

    real du_dy = 0.5 
               * (   (u[Comp::u()][i]  [j+1][k] - u[Comp::u()][i]  [j-1][k]) 
                   + (u[Comp::u()][i+1][j+1][k] - u[Comp::u()][i+1][j-1][k]))
		/(u.dyn(j)+u.dys(j));

    real du_dz = 0.5 
               * (   (u[Comp::u()][i]  [j][k+1] - u[Comp::u()][i]  [j][k-1]) 
                   + (u[Comp::u()][i+1][j][k+1] - u[Comp::u()][i+1][j][k-1]))
		/(u.dzt(k)+u.dzb(k));

    real dv_dx = 0.5 
               * (   (u[Comp::v()][i+1][j]  [k] - u[Comp::v()][i-1][j]  [k]) 
                   + (u[Comp::v()][i+1][j+1][k] - u[Comp::v()][i-1][j+1][k]))
		/(u.dxw(i)+u.dxe(i));

    real dv_dy = (u[Comp::v()][i][j+1][k] - u[Comp::v()][i][j][k]) / u.dyc(j);

    real dv_dz = 0.5 
               * (   (u[Comp::v()][i][j]  [k+1] - u[Comp::v()][i][j]  [k-1]) 
                   + (u[Comp::v()][i][j+1][k+1] - u[Comp::v()][i][j+1][k-1]))
		/(u.dzt(k)+u.dzb(k));

    real dw_dx = 0.5 
               * (   (u[Comp::w()][i+1][j][k]   - u[Comp::w()][i-1][j][k]) 
                   + (u[Comp::w()][i+1][j][k+1] - u[Comp::w()][i-1][j][k+1]))
		/(u.dxw(i)+u.dxe(i));

    real dw_dy = 0.5 
               * (   (u[Comp::w()][i][j+1][k]   - u[Comp::w()][i][j-1][k]  ) 
                + (u[Comp::w()][i][j+1][k+1] - u[Comp::w()][i][j-1][k+1]))
		/(u.dyn(j)+u.dys(j));

    real dw_dz = (u[Comp::w()][i][j][k+1] - u[Comp::w()][i][j][k]) / u.dzc(k);

    /*=====================*/
    real s11 = 0.5*(du_dx + du_dx);
    real s12 = 0.5*(du_dy + dv_dx);
    real s13 = 0.5*(du_dz + dw_dx);
    /*---------------------*/
    real s21 = s12;
    real s22 = 0.5*(dv_dy + dv_dy);
    real s23 = 0.5*(dv_dz + dw_dy);
    /*---------------------*/
    real s31 = s13;
    real s32 = s23;
    real s33 = 0.5*(dw_dz + dw_dz);
    /*=====================*/

    /*=====================*/
    real o11 =  0.0;
    real o12 =  0.5*(du_dy - dv_dx);
    real o13 =  0.5*(du_dz - dw_dx);
    /*---------------------*/
    real o21 = -o12;
    real o22 =  0.0;
    real o23 =  0.5*(dv_dz - dw_dy);
    /*---------------------*/
    real o31 = -o13;
    real o32 = -o23;
    real o33 =  0.0;
    /*=====================*/

    real ss = s11*s11 + s12*s12 + s13*s13 +
              s21*s21 + s22*s22 + s23*s23 +
              s31*s31 + s32*s32 + s33*s33;  

    real oo = o11*o11 + o12*o12 + o13*o13 +
              o21*o21 + o22*o22 + o23*o23 +
              o31*o31 + o32*o32 + o33*o33;  

    (*q)[i][j][k] = 0.5 * ( oo - ss );
  }

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_get_q.cpp,v 1.7 2012/10/08 13:36:48 draksler_m Exp $'/
+-----------------------------------------------------------------------------*/
