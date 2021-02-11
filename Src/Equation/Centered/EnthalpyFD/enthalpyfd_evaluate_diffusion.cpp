#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$
*         or adds diffusion to right hand side.
*******************************************************************************/
void EnthalpyFD::evaluate_diffusion(const Old old, const Scalar * diff_eddy) {

  /* initialize: get time stepping coefficient */
  real tscn = diff_ts.N(); /* 1.0 for fully implicit */
  real tscm =  diff_ts.Nm1(); /* 0.5 for c.n. 0.0 for fully implicit */

  ResistEval re;
  if(old==Old::no) {
    assert( tscn > 0.0 );
    re = ResistEval::no;
    ftif = 0.;
  } else {
    assert( tscm > 0.0 );
    re = ResistEval::yes;
  }

  std::vector<StencilPoint> stencil;
  coef_gen coef_m, coef_p;
  real A_m, A_c, A_p;
  real dSm, dSp;
  int ox(0),oy(0),oz(0);

  for_m(m) {
    if(m==Comp::i()) {
      coef_m = &EnthalpyFD::coef_x_m;
      coef_p = &EnthalpyFD::coef_x_p;
      ox = 1; oy = 0; oz = 0;
    } else if((m==Comp::j())) {
      coef_m = &EnthalpyFD::coef_y_m;
      coef_p = &EnthalpyFD::coef_y_p;
      ox = 0; oy = 1; oz = 0;
    } else {
      coef_m = &EnthalpyFD::coef_z_m;
      coef_p = &EnthalpyFD::coef_z_p;
      ox = 0; oy = 0; oz = 1;
    }
      
    for_ijk(i,j,k) {
      real * Ac = &(A.c [i][j][k]);
      real * Ff = &(ftif[i][j][k]);
      real * Am, * Ap;

      /* A_m, A_c, A_p are dummy variables, only fold is changed in the body */
      if(old==Old::yes) {
        Am = &A_m;
        Ac = &A_c;
        Ap = &A_p;
        Ff = &(fold[i][j][k]);
      } else {
        if(m==Comp::i()) {
          Am = &(A.w[i][j][k]);
          Ap = &(A.e[i][j][k]);
          dSm = phi.dSx(Sign::neg(),i,j,k);
          dSp = phi.dSx(Sign::pos(),i,j,k);
        } else if((m==Comp::j())) {
          Am = &(A.s[i][j][k]);
          Ap = &(A.n[i][j][k]);
          dSm = phi.dSy(Sign::neg(),i,j,k);
          dSp = phi.dSy(Sign::pos(),i,j,k);
        } else {
          Am = &(A.b[i][j][k]);
          Ap = &(A.t[i][j][k]);
          dSm = phi.dSz(Sign::neg(),i,j,k);
          dSp = phi.dSz(Sign::pos(),i,j,k);
        }
      }

      /* center in fluid */
      if(dom->ibody().on(i,j,k)) {
        cell_diffusion_fluid(m,i,j,k,ox,oy,oz,phi.xc(i),
                             coef_m,coef_p,re,old,
                             tscn,tscm,stencil,
                             *Am,*Ac,*Ap,*Ff,
                             diff_eddy);
      } else {
        cell_diffusion_solid(m,i,j,k,ox,oy,oz,phi.xc(i),
                             coef_m,coef_p,dSm,dSp,re,old,
                             tscn,tscm,
                             *Am,*Ac,*Ap,*Ff,
                             diff_eddy);
      } /* solid vs fluid */

    } /* ijk */
  } /* m */

  /*-------------------------------+
  |  non-conductive immersed body  |
  +-------------------------------*/
  if(!solid()&&dom->ibody().nccells() > 0) {
    boil::oout<<"EFD: Conduction without solid. "
              <<"Underdevelopment. Exiting."<<boil::endl;
    exit(0); 

#if 0
    if(old==Old::no) {
      for_ijk(i,j,k)
        if( dom->ibody().off(i,j,k) ) {
          A.c[i][j][k]  = 1.0e-12;
          A.w[i][j][k]  = 0.0;
          A.e[i][j][k]  = 0.0;
          A.s[i][j][k]  = 0.0;
          A.n[i][j][k]  = 0.0;
          A.b[i][j][k]  = 0.0;
          A.t[i][j][k]  = 0.0;
          A.ci[i][j][k] = 1.0;
          fnew[i][j][k] = 0.0;
        }
    }
#endif
  } /* is there an immersed body */

  if(old==Old::no) {
    A.c.exchange();
    A.w.exchange();
    A.e.exchange();
    A.s.exchange();
    A.n.exchange();
    A.b.exchange();
    A.t.exchange();
  }

  return;
}
