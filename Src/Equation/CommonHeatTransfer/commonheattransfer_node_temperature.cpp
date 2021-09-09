#include "commonheattransfer.h"

/******************************************************************************/
void CommonHeatTransfer::calculate_node_temperature(const Scalar * diff_eddy) {
/***************************************************************************//*** 
*  \brief calculate node temperature at the solid/fluid boundary 
*******************************************************************************/
  for_m(m)
    bndtpr_flu(m) = boil::unreal;
  for_m(m)
    bndtpr_sol(m) = boil::unreal;

  Sign dmmy; /* dummy sign */

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<topo->domain()->ibody().nccells(); cc++) {
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    topo->domain()->ibody().ijk(cc,&i,&j,&k);

    Comp m = Comp::i();

    /* west is in wall */
    if(topo->domain()->ibody().off(i-1,j,k)) {

      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i-1,j,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i-1][j][k];

      real len_1 = 0.5*tpr.dxc(i);
      real len_2 = 0.5*tpr.dxc(i-1);

      real res_1 = len_1/lam_1;
      real res_2a = len_2/lam_2;
      real res_2b = wall_resistance(i,j,k);
      real res_2 = res_2a + res_2b;
      real qmult = res_2a / res_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::neg(),m,i,j,k,Old::no)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         len_1 = 
           distance_int_x(Sign::pos(),i-1,j,k,tpr_1,dmmy,ResistEval::no,Old::no)
           -len_2;

         res_1 = len_1/lam_1;
         /* only liquid resistance is considered */
         if(use_int_resist) {
           if(topo->below_interface(i,j,k)) {
             res_1 += int_resistance_liq(i,j,k);
           }
         }

      }

      bndtpr_flu[m][i][j][k] = temperature_node(
                                 dirac_wall_source(i,j,k)*qmult,
                                 res_1, tpr_1,
                                 res_2, tpr_2);
      bndtpr_sol[m][i][j][k] = temperature_node(
                                 dirac_wall_source(i,j,k),
                                 res_2b, bndtpr_flu[m][i][j][k],
                                 res_2a, tpr_2);

    }

    /* east is in wall */
    if(topo->domain()->ibody().off(i+1,j,k)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i+1,j,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i+1][j][k];

      real len_1 = 0.5*tpr.dxc(i);
      real len_2 = 0.5*tpr.dxc(i+1);

      real res_1 = len_1/lam_1;
      real res_2a = len_2/lam_2;
      real res_2b = wall_resistance(i,j,k);
      real res_2 = res_2a + res_2b;
      real qmult = res_2a / res_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::pos(),m,i,j,k,Old::no)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         len_1 = 
           distance_int_x(Sign::neg(),i+1,j,k,tpr_1,dmmy,ResistEval::no,Old::no)
           -len_2;

         res_1 = len_1/lam_1;
         /* only liquid resistance is considered */
         if(use_int_resist) {
           if(topo->below_interface(i,j,k)) {
             res_1 += int_resistance_liq(i,j,k);
           }
         }

      }

      bndtpr_flu[m][i+1][j][k] = temperature_node(
                                   dirac_wall_source(i,j,k)*qmult,
                                   res_1, tpr_1,
                                   res_2, tpr_2);
      bndtpr_sol[m][i+1][j][k] = temperature_node(
                                   dirac_wall_source(i,j,k),
                                   res_2b, bndtpr_flu[m][i+1][j][k],
                                   res_2a, tpr_2);

    }

    m = Comp::j();

    /* south is in wall */
    if(topo->domain()->ibody().off(i,j-1,k)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j-1,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j-1][k];

      real len_1 = 0.5*tpr.dyc(j);
      real len_2 = 0.5*tpr.dyc(j-1);

      real res_1 = len_1/lam_1;
      real res_2a = len_2/lam_2;
      real res_2b = wall_resistance(i,j,k);
      real res_2 = res_2a + res_2b;
      real qmult = res_2a / res_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::neg(),m,i,j,k,Old::no)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         len_1 = 
           distance_int_y(Sign::pos(),i,j-1,k,tpr_1,dmmy,ResistEval::no,Old::no)
           -len_2;

         res_1 = len_1/lam_1;
         /* only liquid resistance is considered */
         if(use_int_resist) {
           if(topo->below_interface(i,j,k)) {
             res_1 += int_resistance_liq(i,j,k);
           }
         }

      }

      bndtpr_flu[m][i][j][k] = temperature_node(
                                 dirac_wall_source(i,j,k)*qmult,
                                 res_1, tpr_1,
                                 res_2, tpr_2);
      bndtpr_sol[m][i][j][k] = temperature_node(
                                 dirac_wall_source(i,j,k),
                                 res_2b, bndtpr_flu[m][i][j][k],
                                 res_2a, tpr_2);

    }

    /* north is in wall */
    if(topo->domain()->ibody().off(i,j+1,k)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j+1,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j+1][k];

      real len_1 = 0.5*tpr.dyc(j);
      real len_2 = 0.5*tpr.dyc(j+1);

      real res_1 = len_1/lam_1;
      real res_2a = len_2/lam_2;
      real res_2b = wall_resistance(i,j,k);
      real res_2 = res_2a + res_2b;
      real qmult = res_2a / res_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::pos(),m,i,j,k,Old::no)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         len_1 = 
           distance_int_y(Sign::neg(),i,j+1,k,tpr_1,dmmy,ResistEval::no,Old::no)
           -len_2;

         res_1 = len_1/lam_1;
         /* only liquid resistance is considered */
         if(use_int_resist) {
           if(topo->below_interface(i,j,k)) {
             res_1 += int_resistance_liq(i,j,k);
           }
         }

      }

      bndtpr_flu[m][i][j+1][k] = temperature_node(
                                   dirac_wall_source(i,j,k)*qmult,
                                   res_1, tpr_1,
                                   res_2, tpr_2);
      bndtpr_sol[m][i][j+1][k] = temperature_node(
                                   dirac_wall_source(i,j,k),
                                   res_2b, bndtpr_flu[m][i][j+1][k],
                                   res_2a, tpr_2);

    }

    m = Comp::k();

    /* bottom is in wall */
    if(topo->domain()->ibody().off(i,j,k-1)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j,k-1,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j][k-1];

      real len_1 = 0.5*tpr.dzc(k);
      real len_2 = 0.5*tpr.dzc(k-1);

      real res_1 = len_1/lam_1;
      real res_2a = len_2/lam_2;
      real res_2b = wall_resistance(i,j,k);
      real res_2 = res_2a + res_2b;
      real qmult = res_2a / res_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::neg(),m,i,j,k,Old::no)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         len_1 = 
           distance_int_z(Sign::pos(),i,j,k-1,tpr_1,dmmy,ResistEval::no,Old::no)
           -len_2;

         res_1 = len_1/lam_1;
         /* only liquid resistance is considered */
         if(use_int_resist) {
           if(topo->below_interface(i,j,k)) {
             res_1 += int_resistance_liq(i,j,k);
           }
         }

      }

      bndtpr_flu[m][i][j][k] = temperature_node(
                                 dirac_wall_source(i,j,k)*qmult,
                                 res_1, tpr_1,
                                 res_2, tpr_2);
      bndtpr_sol[m][i][j][k] = temperature_node(
                                 dirac_wall_source(i,j,k),
                                 res_2b, bndtpr_flu[m][i][j][k],
                                 res_2a, tpr_2);

    }

    /* top is in wall */
    if(topo->domain()->ibody().off(i,j,k+1)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j,k+1,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j][k+1];

      real len_1 = 0.5*tpr.dzc(k);
      real len_2 = 0.5*tpr.dzc(k+1);

      real res_1 = len_1/lam_1;
      real res_2a = len_2/lam_2;
      real res_2b = wall_resistance(i,j,k);
      real res_2 = res_2a + res_2b;
      real qmult = res_2a / res_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::pos(),m,i,j,k,Old::no)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         len_1 = 
           distance_int_z(Sign::neg(),i,j,k+1,tpr_1,dmmy,ResistEval::no,Old::no)
           -len_2;

         res_1 = len_1/lam_1;
         /* only liquid resistance is considered */
         if(use_int_resist) {
           if(topo->below_interface(i,j,k)) {
             res_1 += int_resistance_liq(i,j,k);
           }
         }

      }

      bndtpr_flu[m][i][j][k+1] = temperature_node(
                                   dirac_wall_source(i,j,k)*qmult,
                                   res_1, tpr_1,
                                   res_2, tpr_2);
      bndtpr_sol[m][i][j][k+1] = temperature_node(
                                   dirac_wall_source(i,j,k),
                                   res_2b, bndtpr_flu[m][i][j][k+1],
                                   res_2a, tpr_2);

    }

  }

  return;
}
