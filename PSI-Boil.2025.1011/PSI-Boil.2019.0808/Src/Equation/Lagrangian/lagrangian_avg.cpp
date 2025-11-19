#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::avg() {
/*---------------------------------------+
|  compute volume averaged velocity      |
|  this function I called from main.cpp  |
|  to output the particles' data         |
+---------------------------------------*/

  for_p(p) {

    real u_avg = 0.0;
    real v_avg = 0.0;
    real w_avg = 0.0;

    for_m(m) { 
      real sum = 0.0;
       real den = 0.0;
      real vol = 0.0;
      int ii=0; if(m == Comp::i()) ii++;
      int jj=0; if(m == Comp::j()) jj++;
      int kk=0; if(m == Comp::k()) kk++;

      for_pijk(p, i, j, k) { 
        const real c_a = interface_fraction(i,    j,    k,    p, continuous);
        const real c_b = interface_fraction(i-ii, j-jj, k-kk, p, continuous);
        const real c_av = 0.5*(c_a + c_b);
        const real d_frac = lagrangian_fraction(c_av);
        sum += d_frac * (*u)[m][i][j][k] * (*u).dV(m,i,j,k); 
        vol +=                    d_frac * (*u).dV(m,i,j,k); 
      }
      boil::cart.sum_real(& sum);
      boil::cart.sum_real(& vol);
      if(m==Comp::u()) u_avg = sum / vol;
      if(m==Comp::v()) v_avg = sum / vol;
      if(m==Comp::w()) w_avg = sum / vol;
    }

    boil::oout <<"t_z, "  << p << " " << time->current_time() + time->dt() << " " 
                          << particles[p].z() <<boil::endl;

    boil::oout <<"t_x, "  << p << " " << time->current_time() + time->dt() << " "
                          << particles[p].x() <<boil::endl;

    boil::oout <<"t_y, "  <<  p << " " <<time->current_time() + time->dt() << " " 
                          << particles[p].y() <<boil::endl;

    boil::oout <<"analy_w, " << p << " " << time->current_time() + time->dt()
                             << " " << particles[p].w() <<boil::endl;

    boil::oout <<"avg_w, "   << p << " " << time->current_time() + time->dt()
                             << " "<< w_avg <<boil::endl;

    boil::oout <<"analy_u, " << p << " " << time->current_time() + time->dt() 
                             << " " << particles[p].u() <<boil::endl;

    boil::oout <<"avg_u, "   << p << " " << time->current_time() + time->dt()
                             << " "<< u_avg  <<boil::endl;

    boil::oout <<"analy_v, " << p << " " << time->current_time() + time->dt() 
                             << " " << particles[p].v() <<boil::endl;

    boil::oout <<"avg_v, " << p << " " <<time->current_time() + time->dt()
                           << " " << v_avg <<boil::endl;
  }

  OPR("lagrangian_avg.cpp");
  getchar();

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_avg.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
