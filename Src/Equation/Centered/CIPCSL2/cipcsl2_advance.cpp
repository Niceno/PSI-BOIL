#include "cipcsl2.h"
#include <cmath>
#include <iomanip>
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::advance() {
/***************************************************************************//**
*  \brief Compute convective equation by CIPCSL2 method.
*******************************************************************************/

#ifdef DEBUG
  boil::oout<<"cipcsl2_advance:start\n";
#endif

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  //real max_fext = -1.0e+300;
  for_aijk(i,j,k){   //must be aijk for insert boundary
    sclr[i][j][k]=clr[i][j][k]+time->dt()*fext[i][j][k];
    //if (max_fext<fext[i][j][k]) max_fext=fext[i][j][k];
  }
  //boil::cart.max_real(&max_fext);

  //sclr.bnd_update();
  bdcond(sclr);
  sclr.exchange_all();
  update_node(sclr);

  /*---------------------------------------------------------------------+ 
  | calculate distance function before convection()                      |
  |  in case that phase change takes place                               |
  |  because convection() generate ad-hoc interface due to oscillation.  |
  +---------------------------------------------------------------------*/
  color_minmax();                  // calcumate maxval and minval
#if 0
  if( nredist>=1 ) {               // if sharpen() will be used
    if ( maxval()>0.5 && minval()<0.5 ){  // if interface exists
      if ( max_fext > 0.0 ) {  // if condensation takes place 
        distfunc(clr,24);          // calculate distance function
      }
    }
  }
#endif

  /*-------------------------+
  |  convect color function  |
  +-------------------------*/
  convection();

  /*-------------+
  |  redistance  |
  +-------------*/
  boil::timer.start("cipcsl2 redist");
  color_minmax();
  if( nredist>=1 ) {               // if sharpen() will be used
    if ( maxval()>0.5 && minval()<0.5 ) {  // if interface exists
      //if ( max_fext <= 0.0 ) {  // if no condensation
        distfunc(clr,24);          // calculate distance function
      //}
      if(time->current_step()%nredist==0) redist(localSharpen);
    }
  }
  boil::timer.stop("cipcsl2 redist");

#ifdef IB
  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
  ib_ext_scalar(clr);
  ib_bdcond(clr); 
#endif

  /*------------------+
  |  copy clr to phi  |
  +------------------*/
  for_aijk(i,j,k)
    phi[i][j][k]=clr[i][j][k];

#if 1
  for_ijk(i,j,k){
    if(fabs(phi[i][j][k]-phisurf)<eps_clr) {
      phi[i][j][k]=phisurf+copysign(1.0,phi[i][j][k]-phisurf)*eps_clr;
    }
  }
#endif

#ifdef DEBUG
  boil::oout<<"cipcsl2_advance:end\n";
#endif

#if 0
  plot_f("f.dat");
  plot_sigx("sigx.dat");
  plot_sigy("sigy.dat");
  plot_sigz("sigz.dat");
  plot_sxyz("sxyzi.dat",Comp::i());
  plot_sxyz("sxyzj.dat",Comp::j());
  plot_sxyz("sxyzk.dat",Comp::k());
  boil::plot->plot((*u),clr,"uvw-clr", time->current_step());
  exit(0);
#endif

  return;
}
