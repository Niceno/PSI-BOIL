#include "enthalpyfdadiabatic.h"

/***************************************************************************//**
*  \brief calculate heat flux on wall
*******************************************************************************/
real EnthalpyFDAdiabatic::hflux_wall(const Scalar & val, const Dir din
                          , const Scalar * diff_eddy) {
  //std::cout<<"hflux_wall: "<<din<<"\n";

  Formula form;

  real hflux=0.0; //average heat flux [W/m2]
  real areaw=0.0; //area of wall [m2]

  for( int b=0; b<val.bc().count(); b++ ) {

    if (val.bc().type_decomp(b)) continue;

    /*------------------------+ 
    |  dirichlet (and inlet)  |
    +------------------------*/
    if( val.bc().type(b) == BndType::dirichlet()) {

      Dir d = phi.bc().direction(b);

      if(d != Dir::undefined()) {
        if(d == din) {

          int iof=0, jof=0, kof=0;

          if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
          if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
          if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;
    
          for_vijk(phi.bc().at(b),i,j,k){
            real area = fabs(iof)*val.dSx(i,j,k)
                      + fabs(jof)*val.dSy(i,j,k)
                      + fabs(kof)*val.dSz(i,j,k);
            areaw += area;
            real lc=lambdav;
            real cp_mass = cpv/rhov;
            if((*clr)[i+iof][j+jof][k+kof]>=0.5){
              lc = lambdal;
              cp_mass = cpl/rhol;
            }
            if(diff_eddy){
              lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
            }
            real alen = fabs(val.xc(i)-val.xc(i+iof))
                      + fabs(val.yc(j)-val.yc(j+jof))
                      + fabs(val.zc(k)-val.zc(k+kof));
            hflux += lc*(val[i][j][k]-val[i+iof][j+jof][k+kof])/alen*area;
          }
        }
      }
    }
  }
  boil::cart.sum_real(&hflux);
  boil::cart.sum_real(&areaw);
  //std::cout<<"areaw= "<<areaw<<"\n";
  if(areaw==0){
    return(0.0);
  } else {
    return (hflux/areaw);
  }
}
