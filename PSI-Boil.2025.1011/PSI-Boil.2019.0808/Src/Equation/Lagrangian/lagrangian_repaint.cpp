#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::repaint() {

  boil::timer.start("lagrangian repaint");

  /* set color function to one, i.e. the value of the continuous phase */
  *this = continuous;

  p_id = OFF + boil::milli; /* p_id is the particle id */

  for_p(p) {

    for_pijk(p,i,j,k) {

      val[i][j][k] = interface_fraction(i,j,k,p, val[i][j][k]);
        
      real color = interface_fraction(i,j,k,p,continuous);
      real d_frac = lagrangian_fraction(color);
      if(d_frac >= 0.5) p_id[i][j][k] = p + boil::milli;

    } /* inside particle's box */

  } /* p */

  /* refresh buffers */
  exchange_all();
  p_id.exchange_all();
 
  boil::timer.stop("lagrangian repaint");

}

/******************************************************************************/
void Lagrangian::repaint_box(int p) {

  boil::timer.start("lagrangian repaint_box");

  for_pijk(p,i,j,k) {
    val[i][j][k] = interface_fraction(i,j,k,p, val[i][j][k]);

    real color = interface_fraction(i,j,k,p,continuous);
    real d_frac = lagrangian_fraction(color);
    if(d_frac >= 0.5) p_id[i][j][k] = p + boil::milli;
  } 

  /* refresh buffers */
  exchange_all();
  p_id.exchange_all();

  boil::timer.stop("lagrangian repaint_box");
}

/******************************************************************************/
void Lagrangian::repaint_erase(int p) {

  boil::timer.start("lagrangian repaint_erase");

  for_pijk(p,i,j,k) {

    real color = interface_fraction(i,j,k,p,continuous);
    real d_frac = lagrangian_fraction(color);
    if(d_frac >= 0.5) p_id[i][j][k] = p + boil::milli;

    val[i][j][k] = continuous;     
  } 

  /* refresh buffers */
  exchange_all();
  p_id.exchange_all();

  boil::timer.stop("lagrangian repaint_erase");
}

/******************************************************************************/
void Lagrangian::correct_id(int par) {

  for (int p = par ; p < size(); p++) {

    for_pijk(p,i,j,k) {
      real color = interface_fraction(i,j,k,p,continuous);
      real d_frac = lagrangian_fraction(color);
      if(d_frac >= 0.5) p_id[i][j][k] = p + boil::milli;
    }
  }

  /* refresh buffers */
   p_id.exchange_all();
}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_repaint.cpp,v 1.6 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
