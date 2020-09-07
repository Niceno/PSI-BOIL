#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief calculate heat flux on wall
*******************************************************************************/
real EnthalpyFD::hflux_wall(const Scalar & val, const Dir din
                          , const Scalar * diff_eddy) {
  //std::cout<<"hflux_wall: "<<din<<"\n";

  Formula form;

  real hflux=0.0; //average heat flux [W/m2]
  real areaw=0.0; //area of wall [m2]
  int i,j,k;

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
      	  int of(0);
	  Comp mcomp;
          Sign sig = Sign::undefined();

          if(d == Dir::imin()) { iof++; mcomp = Comp::i(); of = -1; sig = Sign::neg(); }
          if(d == Dir::imax()) { iof--; mcomp = Comp::i(); of = +1; sig = Sign::pos(); }
          if(d == Dir::jmin()) { jof++; mcomp = Comp::j(); of = -1; sig = Sign::neg(); }
          if(d == Dir::jmax()) { jof--; mcomp = Comp::j(); of = +1; sig = Sign::pos(); } 
          if(d == Dir::kmin()) { kof++; mcomp = Comp::k(); of = -1; sig = Sign::neg(); }
          if(d == Dir::kmax()) { kof--; mcomp = Comp::k(); of = +1; sig = Sign::pos(); }
    
          for_vijk(phi.bc().at(b),i,j,k){
            real area = fabs(iof)*val.dSx(sig,i,j,k)
                      + fabs(jof)*val.dSy(sig,i,j,k)
                      + fabs(kof)*val.dSz(sig,i,j,k);
            real alen = fabs(val.xc(i)-val.xc(i+iof))
                      + fabs(val.yc(j)-val.yc(j+jof))
                      + fabs(val.zc(k)-val.zc(k+kof));

            areaw += area;
            if       (dom->ibody().off(i,j,k)) {
              real lc = solid()->lambda(i,j,k);
              hflux += lc*(val[i][j][k]-val[i+iof][j+jof][k+kof])/alen*area;

            } else if(!interface(sig,mcomp,i+iof,j+jof,k+kof)) {
              real lc=lambdav;
              real cp_mass = cpv/rhov;
              if(topo->above_interface(i+of,j+of,k+of)) {
                lc = lambdal;
                cp_mass = cpl/rhol;
              }
              if(diff_eddy){
                lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
              }
              hflux += lc*(val[i][j][k]-val[i+iof][j+jof][k+kof])/alen*area;

	    } else {
              real lc=lambdal;
              real cp_mass = cpl/rhol;
              if(topo->above_interface(i+of,j+of,k+of)) {
                lc = lambdav;
                cp_mass = cpv/rhov;
              }
              if(diff_eddy){
                lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
              }
              real ts;
              if       (mcomp==Comp::i()) {
                alen = fabs(val.xc(i)-val.xc(i+iof))
                     - distance_int_x(sig,i+iof,j+jof,k+kof,ts);
              } else if(mcomp==Comp::j()) {
                alen = fabs(val.yc(j)-val.yc(j+jof))
                     - distance_int_y(sig,i+iof,j+jof,k+kof,ts);
              } else {
                alen = fabs(val.zc(k)-val.zc(k+kof))
                     - distance_int_z(sig,i+iof,j+jof,k+kof,ts);
              }
              hflux += (val[i][j][k]-ts)/(alen/lc+htwallmodel.near_wall_resist)*area;
            } /* interface */
          } /* vijk */
        } /* d == din */
      } /* d != undefined */
    } /* bc == dirichlet */
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

