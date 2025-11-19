#include "floodfill.h"
#include "../../Plot/plot.h"

/***************************************************************************//**
*  computes the properties of each tmp calc ID region    
*  - center of mass of each region x,y,z
*  - volume of each region
*  - average velocity of each region u,v,w
*******************************************************************************/
void Floodfill::get_region_info(int ptmpcalcID, int ntmpcalcID) {

  int size_lookups = ptmpcalcID - ntmpcalcID + 1;
  flookup_cvol.resize(size_lookups);
  flookup_x.resize(size_lookups);
  flookup_y.resize(size_lookups);
  flookup_z.resize(size_lookups);
  flookup_u.resize(size_lookups);
  flookup_v.resize(size_lookups);
  flookup_w.resize(size_lookups);

  /* to ensure reassigned to 0, this is necessary */
  for (int i=0; i<flookup_cvol.size(); i++) {
    flookup_cvol[i] = 0;
    flookup_x[i] = 0.;
    flookup_y[i] = 0.;
    flookup_z[i] = 0.;
    flookup_u[i] = 0.;
    flookup_v[i] = 0.;
    flookup_w[i] = 0.;
  }

  pt_flookup_cvol = flookup_cvol.data();
  pt_flookup_x   = flookup_x.data();
  pt_flookup_y   = flookup_y.data();
  pt_flookup_z   = flookup_z.data();
  pt_flookup_u   = flookup_u.data();
  pt_flookup_v   = flookup_v.data();
  pt_flookup_w   = flookup_w.data();
  pt_flookup_cvol        -= ntmpcalcID;
  pt_flookup_x           -= ntmpcalcID;
  pt_flookup_y           -= ntmpcalcID;
  pt_flookup_z           -= ntmpcalcID;
  pt_flookup_u           -= ntmpcalcID;
  pt_flookup_v           -= ntmpcalcID;
  pt_flookup_w           -= ntmpcalcID;

  /*----------------------------------------------------+
  |  Calculate the volume, center of mass, and avg uvw  |
  |  for each of the final regions                      |
  +----------------------------------------------------*/
  for_vijk(rgnid,i,j,k) {
    int rid = int(rgnid[i][j][k]);
    pt_flookup_cvol[rid]++;
    pt_flookup_x[rid] += rgnid.xc(i);
    pt_flookup_y[rid] += rgnid.yc(j);
    pt_flookup_z[rid] += rgnid.zc(k);
    pt_flookup_u[rid] += 0.5*((*uvw)[Comp::u()][i][j][k]+(*uvw)[Comp::u()][i+1][j][k]);
    pt_flookup_v[rid] += 0.5*((*uvw)[Comp::v()][i][j][k]+(*uvw)[Comp::v()][i][j+1][k]);
    pt_flookup_w[rid] += 0.5*((*uvw)[Comp::w()][i][j][k]+(*uvw)[Comp::w()][i][j][k+1]);
  }
  boil::cart.sum_int_n(&flookup_cvol[0], size_lookups);
  boil::cart.sum_real_n(&flookup_x[0], size_lookups);
  boil::cart.sum_real_n(&flookup_y[0], size_lookups);
  boil::cart.sum_real_n(&flookup_z[0], size_lookups);
  boil::cart.sum_real_n(&flookup_u[0], size_lookups);
  boil::cart.sum_real_n(&flookup_v[0], size_lookups);
  boil::cart.sum_real_n(&flookup_w[0], size_lookups);
  for (int i=ntmpcalcID; i<=ptmpcalcID; i++) {
    if (pt_flookup_cvol[i]) {  
      real volcells = pt_flookup_cvol[i];
      pt_flookup_x[i] /= volcells;
      pt_flookup_y[i] /= volcells;
      pt_flookup_z[i] /= volcells;
      pt_flookup_u[i] /= volcells;
      pt_flookup_v[i] /= volcells;
      pt_flookup_w[i] /= volcells;
    }
  }

}
