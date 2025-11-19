#include "vof.h"

/******************************************************************************/
real VOF::translate_v(const int i, const int j, const int k,
                      const real dx, const real dy, const real dz,
                      const real fx, const real fy, const real fz,
                      real & nnx, real & nny, real & nnz, real & naa,
                      const Scalar & scp) const {
/***************************************************************************//**
*  \brief Calculate volume fraction in an arbitrary cell.
*         i,j,k: original cell
*         dx,...: translation vector (new-old coords)
*         fx,...: scaling vector (dnew/dold)
*         nnx...: new normal vector and alpha
******************************************************************************/

  /* input */
  real scpscp = scp[i][j][k];
  real alphaval = nalpha[i][j][k];

  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -nx[i][j][k];
  real vn2 = -ny[i][j][k];
  real vn3 = -nz[i][j][k];

  /* surface normal is not affected by translation */
  nnx = -vn1;
  nny = -vn2;
  nnz = -vn3;

  /* normalized normal */
  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real denom = vm1+vm2+vm3;

  /* degenerate cases */
  if(!boil::realistic(alphaval)||denom<boil::pico) {
    return real(scpscp>phisurf);
  }

  /* normalized values */
  alphaval /= denom;
  vm1 /= denom;
  vm2 /= denom;
  vm3 /= denom;

  /* new m's */
  real wm1 = vm1*fx;
  real wm2 = vm2*fy;
  real wm3 = vm3*fz;

  real denom2 = wm1+wm2+wm3;

  wm1 /= denom2;
  wm2 /= denom2;
  wm3 /= denom2;

  /* translation effect */
  real d1 = signum(dx/scp.dxc(i),vn1);
  real d2 = signum(dy/scp.dyc(j),vn2);
  real d3 = signum(dz/scp.dzc(k),vn3);

  real translation_sum = d1*vm1+d2*vm2+d3*vm3;

  /* new normalized alpha to be used */
  real alp2 = (alphaval - 0.5 - translation_sum)/denom2 + 0.5;

  /* new unnormalized alpha to be stored */
  naa = alp2*denom2;

  /* volume fraction */
  return calc_v(alp2,wm1,wm2,wm3);

}

/******************************************************************************/
real VOF::extrapolate_v(const int i, const int j, const int k,
                        const int ofx, const int ofy, const int ofz,
                        const real xp, const real yp, const real zp,
                        const Scalar & scp) {
/***************************************************************************//**
*  \brief Simplified version used in update at walls.
*         i,j,k: target cell
*         ii,jj,kk: original cell
*         xp,yp,zp: target cell positions in normalised space
******************************************************************************/

  int ii = i+ofx;
  int jj = j+ofy;
  int kk = k+ofz;
   
  real scpscp = scp[ii][jj][kk];

  /* erroneous interfaces */
  if(scpscp<tol_wall||scpscp-1.0>-tol_wall) {
    return real(scpscp>phisurf);
  } 

  /* unnormalized alpha value */
  real alphaval = nalpha[ii][jj][kk];
          
  /* degenerate case I */
  if(!boil::realistic(alphaval)) {
    return real(scpscp>phisurf);
  } 

  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -nx[ii][jj][kk];
  real vn2 = -ny[ii][jj][kk];
  real vn3 = -nz[ii][jj][kk];

  /* surface normal is not affected by translation */
  nx[i][j][k] = -vn1;
  ny[i][j][k] = -vn2;
  nz[i][j][k] = -vn3;
           
  /* normalized normal */
  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real denom = vm1+vm2+vm3;

  /* degenerate case II */
  if(denom<boil::pico) {
    return real(scpscp>phisurf);
  }
         
  real xpos = xp;
  real ypos = yp;
  real zpos = zp;

  /* mirror boundary cell to the normalized space */
  if(vn1<0) {
    xpos = 1.0-xpos;
  }
  if(vn2<0) {
    ypos = 1.0-ypos;
  }
  if(vn3<0) {
    zpos = 1.0-zpos;
  }

  /* now, we need to translate the coordinate system */
  /* so that x'(center of bnd cell) = (0.5,0.5,0.5)  */
  /* x' = x - x(center of bnd cell) + (0.5,0.5,0.5)  */
  /* this affects alpha value: */
  /* m dot x' = alpha + m dot [(0.5,0.5,0.5) - x(cbc)] = alpha' */
  alphaval += vm1*(0.5-xpos) + vm2*(0.5-ypos) + vm3*(0.5-zpos);

  /* normalized alpha value */
  alphaval /= denom; 
  vm1 /= denom;
  vm2 /= denom;
  vm3 /= denom;

  /* volume fraction */
  return calc_v(alphaval,vm1,vm2,vm3);
}
