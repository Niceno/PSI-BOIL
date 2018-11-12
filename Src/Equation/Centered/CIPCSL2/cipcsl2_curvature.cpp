#include "cipcsl2.h"
#include <iomanip>
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::curvature() {
/***************************************************************************//**
*  \brief Calculate curvature
*******************************************************************************/
  boil::timer.start("cipcsl2 curvature");
#ifdef DEBUG
  std::cout<<"curvature::begin "<<boil::cart.iam()<<"\n";
#endif

  /*-----------+
  |  1st step  |
  +-----------*/
  /* distance function */
  distfunc(clr, 24);

  /* diffusing color function */
  if (!use_dist_for_kappa) {
    for_aijk(i,j,k)
      sclr[i][j][k] = clr [i][j][k];

    sharpen(sclr, 1.5, itsmear, false);
    //smear(sclr);
  }

#ifdef DEBUG
  std::cout<<"tension::1st step end "<<boil::cart.iam()<<"\n";
#endif

  /*-----------+
  |  2nd step  |
  +-----------*/
  /* calculate curvature */
  if (use_dist_for_kappa) {
    curv(dist);
  } else {
    curv(sclr);
  }

#if 0
  if (time->current_step()>=18273 && time->current_step()<=18300) {
  boil::plot->plot(clr,nx,ny,nz, "clr-nx-ny-nz-bef", time->current_step());
  boil::plot->plot(clr,kappa,dist, "clr-kappa-dist-bef", time->current_step());
  }
#endif


  /* wall adhesion */
  //bdcurv(sclr);
  bdcurv(dist);
  //bdcurv(clr);

#if 0
  if (time->current_step()>=18273 && time->current_step()<=18300) {
  boil::plot->plot(clr,nx,ny,nz, "clr-nx-ny-nz-aft", time->current_step());
  boil::plot->plot(clr,kappa,dist, "clr-kappa-dist-aft", time->current_step());
  //exit(0);
  }
#endif

  boil::timer.stop("cipcsl2 curvature");
}
