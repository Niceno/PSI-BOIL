#include "twolevelvector.h"

/* Divergence-free interpolation of velocities (2D, X-Z)  */
/* input: coarse, output: fine, throughput: cs (coarse scalar)
                                            mdot (fine mass source)
                                            fluid (matter)

     ---------     b = boil::BW-1
     | x | y |     coarse coords: i,j,k ; i^ = i-b, k^ = k-b
     ---------     fine coords:
     | q | z |         y - 2*i^   +b, j ,2*k^   +b
     ---------         x - 2*i^-1 +b, j ,2*k^   +b
                       z - 2*i^   +b, j ,2*k^-1 +b
                       q - 2*i^-1 +b, j ,2*k^-1 +b

         T-W     T-E         
      -----------------
    W |       |       |  E
    - |   x   A   y   |  -
    T |       |       |  T
      ----D-------B----
    W |       |       |  E
    - |   q   C   z   |  -
    B |       |       |  B
      -----------------
         B-W     B-E         

    unknown velocities: A, B, C, D
    auxilliary unknown scalar values: x, y, q, z

*/

void TwoLevelVector::divergence_free_interpolate_XZ(
                   const Scalar & cs, const Scalar & mdot,
                   const Matter & fluid) {

  real mult(0.);
  if(fluid.mixture()) {
    const real rhov = fluid.rho(0);
    const real rhol = fluid.rho(1);
    mult = 1./rhov-1./rhol;
  }
  
  for_vijk(cs,i,j,k) {
   
    /* step 1: establish coordinate correspondence */
    const int b = boil::BW-1;
    const int ihat = i-b;
    const int khat = k-b;

    const int i_y = 2*ihat + b;
    const int j_y = j;
    const int k_y = 2*khat + b;

    const int i_x = 2*ihat-1 + b;
    const int j_x = j;
    const int k_x = 2*khat + b;

    const int i_z = 2*ihat + b;
    const int j_z = j;
    const int k_z = 2*khat-1 + b;

    const int i_q = 2*ihat-1 + b;
    const int j_q = j;
    const int k_q = 2*khat-1 + b;

    /* step 2: calculate fine side-velocities */
    /* --- currently: constant velocity field assumed --- */
    Comp m = Comp::u();

    const real u_wt = coarse[m][i  ][j][k];
    const real u_wb = coarse[m][i  ][j][k];

    const real u_et = coarse[m][i+1][j][k];
    const real u_eb = coarse[m][i+1][j][k];

    m = Comp::k();

    const real u_bw = coarse[m][i][j][k  ];
    const real u_be = coarse[m][i][j][k  ];

    const real u_tw = coarse[m][i][j][k+1];
    const real u_te = coarse[m][i][j][k+1];

    /* step 3: calculate flow areas */
    const real S_wt = mdot.dSx(Sign::neg(),i_x,j_x,k_x); 
    const real S_wb = mdot.dSx(Sign::neg(),i_q,j_q,k_q); 

    const real S_a = mdot.dSx(Sign::pos(),i_x,j_x,k_x);
    const real S_c = mdot.dSx(Sign::pos(),i_q,j_q,k_q);

    const real S_et = mdot.dSx(Sign::pos(),i_y,j_y,k_y); 
    const real S_eb = mdot.dSx(Sign::pos(),i_z,j_z,k_z); 

    const real S_bw = mdot.dSz(Sign::neg(),i_q,j_q,k_q); 
    const real S_be = mdot.dSz(Sign::neg(),i_z,j_z,k_z); 

    const real S_d = mdot.dSz(Sign::pos(),i_q,j_q,k_q); 
    const real S_b = mdot.dSz(Sign::pos(),i_z,j_z,k_z); 

    const real S_tw = mdot.dSz(Sign::pos(),i_x,j_x,k_x); 
    const real S_te = mdot.dSz(Sign::pos(),i_y,j_y,k_y); 

    /* step 4: calculate side flow rates */
    const real Q_wt = S_wt*u_wt;
    const real Q_wb = S_wb*u_wb;

    const real Q_et = S_et*u_et;
    const real Q_eb = S_eb*u_eb;

    const real Q_bw = S_bw*u_bw;
    const real Q_be = S_be*u_be;

    const real Q_tw = S_tw*u_tw;
    const real Q_te = S_te*u_te;

    /* step 5: internal flow rates are estimated through averaging */
    const real Q_a = 0.5*(Q_wt+Q_et);
    const real Q_c = 0.5*(Q_wb+Q_eb);

    const real Q_d = 0.5*(Q_bw+Q_tw);
    const real Q_b = 0.5*(Q_be+Q_te);

    /* step 6: calculate matrix coefficients */
    const real C_a = S_a/mdot.dxw(i_y);
    const real C_c = S_c/mdot.dxw(i_z);

    const real C_d = S_d/mdot.dzb(k_x);
    const real C_b = S_b/mdot.dzb(k_y);

    /* step 7: calculate source terms */
    const real F_x = mult*mdot[i_x][j_x][k_x]*mdot.dV(i_x,j_x,k_x);
    const real F_y = mult*mdot[i_y][j_y][k_y]*mdot.dV(i_y,j_y,k_y);
    const real F_z = mult*mdot[i_z][j_z][k_z]*mdot.dV(i_z,j_z,k_z);
    const real F_q = mult*mdot[i_q][j_q][k_q]*mdot.dV(i_q,j_q,k_q);

    /* step 8: solve for pressures x,y,z,q */
    real x, y, z, q;

    interpolate_pressure_solve_2D(x,y,z,q,
                                  F_x,F_y,F_z,F_q,
                                  C_a,C_b,C_c,C_d,
                                  Q_a,Q_b,Q_c,Q_d,
                                  Q_wb,Q_wt,Q_eb,Q_et,
                                  Q_bw,Q_be,Q_tw,Q_te);

    /* step 9: project onto internal fine velocities */
    const real u_a = (Q_a-C_a*(y-x))/S_a;
    const real u_c = (Q_c-C_c*(z-q))/S_c;

    const real u_d = (Q_d-C_d*(x-q))/S_d;
    const real u_b = (Q_b-C_b*(y-z))/S_b;

    /* step 10: assign values to fine vector field */
    m = Comp::u();

    fine[m][i_x  ][j_x][k_x] = u_wt;
    fine[m][i_y  ][j_y][k_y] = u_a;
    fine[m][i_y+1][j_y][k_y] = u_et;

    fine[m][i_q  ][j_q][k_q] = u_wb;
    fine[m][i_z  ][j_z][k_z] = u_c;
    fine[m][i_z+1][j_z][k_z] = u_eb;

    m = Comp::w();

    fine[m][i_q][j_q][k_q  ] = u_bw;
    fine[m][i_x][j_x][k_x  ] = u_d;
    fine[m][i_x][j_x][k_x+1] = u_tw;

    fine[m][i_z][j_z][k_z  ] = u_be;
    fine[m][i_y][j_y][k_y  ] = u_b;
    fine[m][i_y][j_y][k_y+1] = u_te;

  } /* ijk */
}

/* the solution to the poisson equation was obtained with python sympy */
void TwoLevelVector::interpolate_pressure_solve_2D(
                   real & x, real & y, real & z, real & q,
                   const real & F_x, const real & F_y,
                   const real & F_z, const real & F_q,
                   const real & c_a, const real & c_b,
                   const real & c_c, const real & c_d,
                   const real & Q_a, const real & Q_b,
                   const real & Q_c, const real & Q_d,
                   const real & Q_wb, const real & Q_wt,
                   const real & Q_eb, const real & Q_et,
                   const real & Q_bw, const real & Q_be,
                   const real & Q_tw, const real & Q_te) {

  /* three equations with three unknowns determined up to a constant;
     thus, q = 0 */
  x = (  F_x*c_a*c_b + F_x*c_a*c_c + F_x*c_b*c_c + F_y*c_a*c_b
       + F_y*c_a*c_c + F_z*c_a*c_b - Q_a*c_b*c_c + Q_b*c_a*c_c
       + Q_be*c_a*c_b + Q_c*c_a*c_b + Q_d*c_a*c_b + Q_d*c_a*c_c
       + Q_d*c_b*c_c - Q_eb*c_a*c_b - Q_et*c_a*c_b - Q_et*c_a*c_c
       - Q_te*c_a*c_b - Q_te*c_a*c_c - Q_tw*c_a*c_b - Q_tw*c_a*c_c
       - Q_tw*c_b*c_c + Q_wt*c_a*c_b + Q_wt*c_a*c_c + Q_wt*c_b*c_c)
      /(c_a*c_b*c_c + c_a*c_b*c_d + c_a*c_c*c_d + c_b*c_c*c_d);

  y = (  c_b*(c_a + c_d)
       * ( c_b*(c_a*(F_x - Q_a + Q_d - Q_tw + Q_wt)
         + (c_a + c_d)*(F_y + Q_a + Q_b - Q_et - Q_te))
         - (c_a*c_a - (c_a + c_b)*(c_a + c_d))*(F_z - Q_b + Q_be + Q_c - Q_eb))
       - ( c_a*(F_x - Q_a + Q_d - Q_tw + Q_wt)
         + (c_a + c_d)*(F_y + Q_a + Q_b - Q_et - Q_te))
       * (c_b*c_b*(c_a + c_d) + (c_a*c_a - (c_a + c_b)*(c_a + c_d))*(c_b + c_c))
      )
      /(( c_a*c_a - (c_a + c_b)*(c_a + c_d))
        * (c_b*c_b*(c_a + c_d) + (c_a*c_a - (c_a + c_b)*(c_a + c_d))*(c_b + c_c))
      );

  z= (  -c_b
      * ( c_a*(F_x - Q_a + Q_d - Q_tw + Q_wt)
        + (c_a + c_d)*(F_y + Q_a + Q_b - Q_et - Q_te))
      + (c_a*c_a - (c_a + c_b)*(c_a + c_d))*(F_z - Q_b + Q_be + Q_c - Q_eb)
     )
    /(c_b*c_b*(c_a + c_d) + (c_a*c_a - (c_a + c_b)*(c_a + c_d))*(c_b + c_c));

  q = 0.;

  return;
}
