#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::calculate_node_temperature(const Scalar * diff_eddy) {
/***************************************************************************//*** 
*  \brief calculate node temperature at the solid/fluid boundary 
*******************************************************************************/
  for_m(m)
    bndtpr(m) = boil::unreal;

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    Comp m = Comp::i();

    /* west is in wall */
    if(dom->ibody().off(i-1,j,k)) {

      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i-1,j,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i-1][j][k];

      real len_1 = 0.5*phi.dxc(i);
      real len_2 = 0.5*phi.dxc(i-1);

      real res_1 = len_1/lam_1;
      real res_2 = len_2/lam_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::neg(),m,i,j,k)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         //len_1 -= distance_int_x(Sign::neg(),i,j,k,tpr_1);
         len_1 = distance_int_x(Sign::pos(),i-1,j,k,tpr_1)-len_2;
         res_1 = len_1/lam_1 + cht.heat_transfer_wall_model().near_wall_resist;
      }

      bndtpr[m][i][j][k] = cht.temperature_node(
                                cht.heat_transfer_wall_model().dirac_wall_source,
                                res_1, tpr_1,
                                res_2, tpr_2);

    }

    /* east is in wall */
    if(dom->ibody().off(i+1,j,k)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i+1,j,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i+1][j][k];

      real len_1 = 0.5*phi.dxc(i);
      real len_2 = 0.5*phi.dxc(i+1);

      real res_1 = len_1/lam_1;
      real res_2 = len_2/lam_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::pos(),m,i,j,k)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         //len_1 -= distance_int_x(Sign::pos(),i,j,k,tpr_1);
         len_1 = distance_int_x(Sign::neg(),i+1,j,k,tpr_1)-len_2;
         res_1 = len_1/lam_1 + cht.heat_transfer_wall_model().near_wall_resist;
      }

      bndtpr[m][i+1][j][k] = cht.temperature_node(
                                  cht.heat_transfer_wall_model().dirac_wall_source,
                                  res_1, tpr_1,
                                  res_2, tpr_2);

    }

    m = Comp::j();

    /* south is in wall */
    if(dom->ibody().off(i,j-1,k)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j-1,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j-1][k];

      real len_1 = 0.5*phi.dyc(j);
      real len_2 = 0.5*phi.dyc(j-1);

      real res_1 = len_1/lam_1;
      real res_2 = len_2/lam_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::neg(),m,i,j,k)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         //len_1 -= distance_int_y(Sign::neg(),i,j,k,tpr_1);
         len_1 = distance_int_y(Sign::pos(),i,j-1,k,tpr_1)-len_2;
         res_1 = len_1/lam_1 + cht.heat_transfer_wall_model().near_wall_resist;
      }

      bndtpr[m][i][j][k] = cht.temperature_node(
                                cht.heat_transfer_wall_model().dirac_wall_source,
                                res_1, tpr_1,
                                res_2, tpr_2);

    }

    /* north is in wall */
    if(dom->ibody().off(i,j+1,k)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j+1,k,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j+1][k];

      real len_1 = 0.5*phi.dyc(j);
      real len_2 = 0.5*phi.dyc(j+1);

      real res_1 = len_1/lam_1;
      real res_2 = len_2/lam_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::pos(),m,i,j,k)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         //len_1 -= distance_int_y(Sign::pos(),i,j,k,tpr_1);
         len_1 = distance_int_y(Sign::neg(),i,j+1,k,tpr_1)-len_2;
         res_1 = len_1/lam_1 + cht.heat_transfer_wall_model().near_wall_resist;
      }

      bndtpr[m][i][j+1][k] = cht.temperature_node(
                                  cht.heat_transfer_wall_model().dirac_wall_source,
                                  res_1, tpr_1,
                                  res_2, tpr_2);

    }

    m = Comp::k();

    /* bottom is in wall */
    if(dom->ibody().off(i,j,k-1)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j,k-1,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j][k-1];

      real len_1 = 0.5*phi.dzc(k);
      real len_2 = 0.5*phi.dzc(k-1);

      real res_1 = len_1/lam_1;
      real res_2 = len_2/lam_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::neg(),m,i,j,k)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         //len_1 -= distance_int_z(Sign::neg(),i,j,k,tpr_1);
         len_1 = distance_int_z(Sign::pos(),i,j,k-1,tpr_1)-len_2;
         res_1 = len_1/lam_1 + cht.heat_transfer_wall_model().near_wall_resist;
      }

      bndtpr[m][i][j][k] = cht.temperature_node(
                                cht.heat_transfer_wall_model().dirac_wall_source,
                                res_1, tpr_1,
                                res_2, tpr_2);

    }

    /* top is in wall */
    if(dom->ibody().off(i,j,k+1)) {
      real lam_1 = lambda(i,j,k,diff_eddy);
      real lam_2 = lambda(i,j,k+1,diff_eddy);

      real tpr_1 = tpr[i][j][k];
      real tpr_2 = tpr[i][j][k+1];

      real len_1 = 0.5*phi.dzc(k);
      real len_2 = 0.5*phi.dzc(k+1);

      real res_1 = len_1/lam_1;
      real res_2 = len_2/lam_2;

      /* inversion of fluid due to existence of interface */
      if(interface(Sign::pos(),m,i,j,k)) {
         lam_1 = lambda_inv(i,j,k,diff_eddy);
         //len_1 -= distance_int_z(Sign::pos(),i,j,k,tpr_1);
         len_1 = distance_int_z(Sign::neg(),i,j,k+1,tpr_1)-len_2;
         res_1 = len_1/lam_1 + cht.heat_transfer_wall_model().near_wall_resist;
      }

      bndtpr[m][i][j][k+1] = cht.temperature_node(
                                  cht.heat_transfer_wall_model().dirac_wall_source,
                                  res_1, tpr_1,
                                  res_2, tpr_2);

    }

  }

  return;
}
