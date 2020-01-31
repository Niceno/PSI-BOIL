#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::near_wall_model(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief For stability reasons, it is probably better to calculate the phase 
*         change rate from the heat flux in the underlying solid. This avoids
*         a sharp spike in heat flux when the interface crosses the cell cent-
*         re as this is the moment the gradt method normally switches from 1st
*         order to second order. 
*
*         Current implementation only works if the condition is boiling (i.e.
*         the heat transfer path is solid->liquid->interface or solid->vapour)
*         and the heat transfer surface is located in the negative-z direction.
*         This is, of course, easily extendable.
*          
*         This function is rather duplicit with respect to correct_gradt_at_ib
*         and somehow questions the need for the involved gradient extrapolation
*         method. The argument is that the previous functions are generic; this
*         one simplifies the problem using a priori knowledge (e.g. lateral HF
*         is assumed to be negligible).
*******************************************************************************/
  
  /* only if heat transfer occurs in the solid */
  if(solid()) {
    for_ijk(i,j,k) {

      /* z direction */
      Comp m = Comp::w();

      /* bottom is in wall and this is an interfacial cell */
      if(dom->ibody().off(i,j,k-1) && Interface(i,j,k)) {
        real lv = lambdav;
        real ll = lambdal;
        
        if (diff_eddy) {
          lv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
          ll += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
        }

        real tmp_w = bndtpr[m][i][j][k];
        real clrc = clr[i][j][k];
        real len_s = phi.dzb(k) - 0.5*phi.dzc(k);
        real lam_s = solid()->lambda(i,j,k-1);
        real tmp_s = tpr[i][j][k-1];

        /* interface exists between the solid and cell centre */
        if(Interface(-1,m,i,j,k)) {
          real lam_f;
          if(clrc>=clrsurf) { /* note the inversion! */
            /* this condition should never occur!!! */
            lam_f=lv;
          } else {
            lam_f=ll;
          }
          real tmp_f;
          real len_f = distance_z(i,j,k,-1,tmp_f);
          len_f = 0.5*phi.dzc(k) - len_f;

          tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
        }
        assert(boil::realistic(tmp_w));

        /* we assume only liquid plays a role and is 'between' solid and gas */
        tnl[i][j][k] = lam_s/ll*(tmp_w-tmp_s)/len_s * (-1);
        tnv[i][j][k] = 0.0;
      }
    } /* ijk */
  } /* solid */

  return;
}
