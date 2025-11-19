#include "vof.h"

/******************************************************************************/
void VOF::ev_discretize(const Matter * fluid, const ScalarInt & pflag,
                        Matrix & A) {
/***************************************************************************//**
*  \brief Construct the system matrix for the Poisson problem.
*  pflag = -1,1 : solved and inside the solution domain
*  pflag = -2 : solved, at least one face is zero-Dirichlet boundary
*  pflag = +2 : solved, at least one face is zero-Neumann boundary
*  pflag = -3,4: bulk cell forming boundary with -2,+2 cells
*  pflag = -1000: solid cell, always zero-Neumann boundary
*******************************************************************************/

  Comp m;
  real a_w, a_e, a_s, a_n, a_b, a_t;
  real rhom, rhop;

  /* coefficients in i direction (w and e) */
  m = Comp::u();
  for_ijk(i,j,k) {
    A.c[i][j][k] = 0.;
    A.w[i][j][k] = 0.;
    A.e[i][j][k] = 0.;
    if(abs(pflag[i][j][k])<3) {
      if(pflag[i-1][j][k]!=-1000) { /* !wall */
        rhom = fluid->rho(m,i,  j,k);
        a_w = dSx(Sign::neg(),i,j,k);
        A.w[i][j][k] = a_w / dxw(i) / rhom;
      }
      if(pflag[i+1][j][k]!=-1000) { /* !wall */
        rhop = fluid->rho(m,i+1,j,k);
        a_e = dSx(Sign::pos(),i,j,k);
        A.e[i][j][k] = a_e / dxe(i) / rhop;
      }
      A.c[i][j][k] += A.w[i][j][k] + A.e[i][j][k];

      /* Neumann boundary */
      if(pflag[i][j][k]==2) {
        if(pflag[i-1][j][k]==4) {
          A.c[i][j][k] -= A.w[i][j][k];
          A.w[i][j][k] = 0.;
        }
        if(pflag[i+1][j][k]==4) {
          A.c[i][j][k] -= A.e[i][j][k];
          A.e[i][j][k] = 0.;
        }
      }
      /* zero-Dirichlet boundary */
      if(pflag[i][j][k]==-2) {
        if(pflag[i-1][j][k]==-3) {
          A.w[i][j][k] = 0.;
        }
        if(pflag[i+1][j][k]==-3) {
          A.e[i][j][k] = 0.;
        }
      }

    } /* pflag < 3 */
  } /* ijk */

  /* coefficients in j direction (s and n) */
  m = Comp::v();
  for_ijk(i,j,k) {
    A.s[i][j][k] = 0.;
    A.n[i][j][k] = 0.;
    if(abs(pflag[i][j][k])<3) {
      if(pflag[i][j-1][k]!=-1000) { /* !wall */
        rhom = fluid->rho(m,i,j,  k);
        a_s = dSy(Sign::neg(),i,j,k);
        A.s[i][j][k] = a_s / dys(j) / rhom;
      }
      if(pflag[i][j+1][k]!=-1000) { /* !wall */
        rhop = fluid->rho(m,i,j+1,k);
        a_n = dSy(Sign::pos(),i,j,k);
        A.n[i][j][k] = a_n / dyn(j) / rhop;
      }
      A.c[i][j][k] += A.s[i][j][k] + A.n[i][j][k];

      /* Neumann boundary */
      if(pflag[i][j][k]==2) {
        if(pflag[i][j-1][k]==4) {
          A.c[i][j][k] -= A.s[i][j][k];
          A.s[i][j][k] = 0.;
        }
        if(pflag[i][j+1][k]==4) {
          A.c[i][j][k] -= A.n[i][j][k];
          A.n[i][j][k] = 0.;
        }
      }
      /* zero-Dirichlet boundary */
      if(pflag[i][j][k]==-2) {
        if(pflag[i][j-1][k]==-3) {
          A.s[i][j][k] = 0.;
        }
        if(pflag[i][j+1][k]==-3) {
          A.n[i][j][k] = 0.;
        }
      }

    } /* pflag < 3 */
  } /* ijk */

  /* coefficients in k direction (b and t) */
  m = Comp::w();
  for_ijk(i,j,k) {
    A.b[i][j][k] = 0.;
    A.t[i][j][k] = 0.;
    if(abs(pflag[i][j][k])<3) {
      if(pflag[i][j][k-1]!=-1000) { /* !wall */
        rhom = fluid->rho(m,i,j,k);
        a_b = dSz(Sign::neg(),i,j,k);
        A.b[i][j][k] = a_b / dzb(k) / rhom;
      }
      if(pflag[i][j][k+1]!=-1000) { /* !wall */
        rhop = fluid->rho(m,i,j,k+1);
        a_t = dSz(Sign::pos(),i,j,k);
        A.t[i][j][k] = a_t / dzt(k) / rhop;
      }
      A.c[i][j][k] += A.b[i][j][k] + A.t[i][j][k];

      /* Neumann boundary */
      if(pflag[i][j][k]==2) {
        if(pflag[i][j][k-1]==4) {
          A.c[i][j][k] -= A.b[i][j][k];
          A.b[i][j][k] = 0.;
        }
        if(pflag[i][j][k+1]==4) {
          A.c[i][j][k] -= A.t[i][j][k];
          A.t[i][j][k] = 0.;
        }
      }
      /* zero-Dirichlet boundary */
      if(pflag[i][j][k]==-2) {
        if(pflag[i][j][k-1]==-3) {
          A.b[i][j][k] = 0.;
        }
        if(pflag[i][j][k+1]==-3) {
          A.t[i][j][k] = 0.;
        }
      }
      
    } /* pflag < 3 */
  } /* ijk */

  /* inverse of the central coefficient */
  for_ijk(i,j,k) {
    if(abs(pflag[i][j][k])<3) {
      A.c[i][j][k] = -A.c[i][j][k];
      A.ci[i][j][k] = 1.0 / A.c[i][j][k];
    } else {
      A.ci[i][j][k] = 0.;
    }
  }

  return;
}

