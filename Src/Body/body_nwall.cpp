#include "body.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"
#include "../Domain/domain.h"
#include <iomanip>
//#define DEBUG

/******************************************************************************/
void Body::nwall(const Domain & dom) {
/***************************************************************************//**
*  /brief calculate normal vector of wall.
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"body_nwall::start "<<boil::cart.iam()<<"\n";
#endif

  /* define scalar */
  ux = new Scalar(dom); (* ux) = 0.0; (* ux)=bdist->shape();
  uy = new Scalar(dom); (* uy) = 0.0; (* uy)=bdist->shape();
  uz = new Scalar(dom); (* uz) = 0.0; (* uz)=bdist->shape();

  for_vijk((*bdist),i,j,k){
    real nx, ny, nz;
    nx=((*bdist)[i+1][j][k]-(*bdist)[i-1][j][k])/(bdist->dxw(i)+bdist->dxe(i));
    ny=((*bdist)[i][j+1][k]-(*bdist)[i][j-1][k])/(bdist->dys(j)+bdist->dyn(j));
    nz=((*bdist)[i][j][k+1]-(*bdist)[i][j][k-1])/(bdist->dzb(k)+bdist->dzt(k));
    normalize(nx, ny, nz);
    (*ux)[i][j][k] = -nx;
    (*uy)[i][j][k] = -ny;
    (*uz)[i][j][k] = -nz;
  }

  (*ux).exchange_all();
  (*uy).exchange_all();
  (*uz).exchange_all();

#ifdef DEBUG
  std::cout<<"body_nwall::end "<<boil::cart.iam()<<"\n";
#endif

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: body_nwall.cpp,v 1.1 2014/02/03 14:12:34 sato Exp $'/
+-----------------------------------------------------------------------------*/
