#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"
#include <iomanip>
#include <iostream>

/******************************************************************************/
void Domain::statistics(Body * body) {
/*--------------------------------------------------------------+
|  Computes grid statistics.                                    |
+--------------------------------------------------------------*/

  /* browse in "i", "j" and "k" direction */
  min_dx = FLT_MAX;
  max_dx = 0.0;
  min_dy = FLT_MAX;
  max_dy = 0.0;
  min_dz = FLT_MAX;
  max_dz = 0.0;
  min_dV = FLT_MAX;
  min_dV = FLT_MAX;
  max_dV = 0.0;         
  max_ar = 0.0;

  int ibeg(boil::BW), iend(ni()-boil::BW);
  int jbeg(boil::BW), jend(nj()-boil::BW);
  int kbeg(boil::BW), kend(nk()-boil::BW);
  
  for(int i=ibeg; i<iend; i++) 
    for(int j=jbeg; j<jend; j++) 
      for(int k=kbeg; k<kend; k++) {         
        if(body){
          if((*body).off(i,j,k)) continue;
        }
        min_dx = std::min(dxc(i), min_dx);
        max_dx = std::max(dxc(i), max_dx);
        min_dy = std::min(dyc(j), min_dy);
        max_dy = std::max(dyc(j), max_dy);
        min_dz = std::min(dzc(k), min_dz);
        max_dz = std::max(dzc(k), max_dz);
        min_dV = std::min(dV(i,j,k), min_dV);
        max_dV = std::max(dV(i,j,k), max_dV);

        real ratio_xy = std::max(dxc(i)/dyc(j), dyc(j)/dxc(i));
        real ratio_xz = std::max(dxc(i)/dzc(k), dzc(k)/dxc(i));
        real ratio_yz = std::max(dyc(j)/dzc(k), dzc(k)/dyc(j));
        real ratio = boil::maxr(ratio_xy, ratio_xz, ratio_yz);

        max_ar = std::max(max_ar, ratio);
      }

  /* find global exreme values */
  boil::cart.min_real(& min_dx);  boil::cart.max_real(& max_dx);
  boil::cart.min_real(& min_dy);  boil::cart.max_real(& max_dy);
  boil::cart.min_real(& min_dz);  boil::cart.max_real(& max_dz);
  boil::cart.min_real(& min_dV);  boil::cart.max_real(& max_dV);
  boil::cart.max_real(& max_ar);

  min_dxyz = boil::minr(min_dx, min_dy, min_dz);
  max_dxyz = boil::maxr(max_dx, max_dy, max_dz);

  const real min_x = global_min_x();
  const real max_x = global_max_x();
  const real lx    = max_x - min_x;
  const real min_y = global_min_y();
  const real max_y = global_max_y();
  const real ly    = max_y - min_y;
  const real min_z = global_min_z();
  const real max_z = global_max_z();
  const real lz    = max_z - min_z;

  int gnx(gi()-2*boil::BW);
  int gny(gj()-2*boil::BW);
  int gnz(gk()-2*boil::BW);

  char b[64];
  boil::oout << "+=========================";
  boil::oout << "==========================+" << boil::endl;
  boil::oout << "| Grid statistics:        ";
  if(body){
  boil::oout << "  (in fluid domain)       |" << boil::endl;
  } else {
  boil::oout << "                          |" << boil::endl;
  }
  boil::oout << "+-------------------------";
  boil::oout << "--------------------------+" << boil::endl;
  sprintf(b, "| Resolution: nx, ny, nz = %6d, %6d, %6d.  |", 
             gnx, gny, gnz);
  boil::oout << b << boil::endl;
  boil::oout << "+-------------------------";
  boil::oout << "+-------------------------+" << boil::endl;
  sprintf(b, "| min x    = %12.5e ",  min_x ); boil::oout << b;
  sprintf(b, "| min dx   = %12.5e |", min_dx); boil::oout << b << boil::endl;
  sprintf(b, "| max x    = %12.5e ",  max_x ); boil::oout << b;
  sprintf(b, "| max dx   = %12.5e |", max_dx); boil::oout << b << boil::endl;
  sprintf(b, "| min y    = %12.5e ",  min_y ); boil::oout << b;
  sprintf(b, "| min dy   = %12.5e |", min_dy); boil::oout << b << boil::endl;
  sprintf(b, "| max y    = %12.5e ",  max_y ); boil::oout << b;
  sprintf(b, "| max dy   = %12.5e |", max_dy); boil::oout << b << boil::endl;
  sprintf(b, "| min z    = %12.5e ",  min_z ); boil::oout << b;
  sprintf(b, "| min dz   = %12.5e |", min_dz); boil::oout << b << boil::endl;
  sprintf(b, "| max z    = %12.5e ",  max_z ); boil::oout << b;
  sprintf(b, "| max dz   = %12.5e |", max_dz); boil::oout << b << boil::endl;
  boil::oout << "+-------------------------";
  boil::oout << "+-------------------------+" << boil::endl;
  sprintf(b, "| lx       = %12.5e ",     lx ); boil::oout << b;
  sprintf(b, "| min dV   = %12.5e |", min_dV); boil::oout << b << boil::endl;
  sprintf(b, "| ly       = %12.5e ",     ly ); boil::oout << b;
  sprintf(b, "| max dV   = %12.5e |", max_dV); boil::oout << b << boil::endl;
  sprintf(b, "| lz       = %12.5e ",     lz ); boil::oout << b;
  sprintf(b, "| max a.r. = %12.5e |", max_ar); boil::oout << b << boil::endl;
  boil::oout << "+-------------------------";
  boil::oout << "+-------------------------+" << boil::endl;
}
