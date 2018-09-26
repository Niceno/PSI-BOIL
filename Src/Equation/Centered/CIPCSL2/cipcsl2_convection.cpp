#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::convection() {
/***************************************************************************//**
*  \brief Compute convective equation by CIPCSL2 method.
*******************************************************************************/

  boil::timer.start("cipcsl2 convection");
  /* reset sum_outlet */
  sum_outlet=0.0;
  sum_outletm=0.0;

  /* set loop range for convert & invert of clr */
  int ist=si(); int ied=ei();
  int jst=sj(); int jed=ej();
  int kst=sk(); int ked=ek();
  if( clr.bc().type_here(Dir::imin(),BndType::periodic())
    ||dom->coord(Comp::i()) !=0)                      ist=si()-1;
  if( clr.bc().type_here(Dir::imax(),BndType::periodic())
    ||dom->coord(Comp::i()) != dom->dim(Comp::i())-1) ied=ei()+1;
  if( clr.bc().type_here(Dir::jmin(),BndType::periodic())
    ||dom->coord(Comp::j()) !=0)                      jst=sj()-1;
  if( clr.bc().type_here(Dir::jmax(),BndType::periodic())
    ||dom->coord(Comp::j()) != dom->dim(Comp::j())-1) jed=ej()+1;
  if( clr.bc().type_here(Dir::kmin(),BndType::periodic())
    ||dom->coord(Comp::k()) !=0)                      kst=sk()-1;
  if( clr.bc().type_here(Dir::kmax(),BndType::periodic())
    ||dom->coord(Comp::k()) != dom->dim(Comp::k())-1) ked=ek()+1;

  /* convert color-function(clr) to density rho */
  for(int i=ist; i<=ied; i++){
    for(int j=jst; j<=jed; j++){
      for(int k=kst; k<=ked; k++){
        clr[i][j][k]=clr[i][j][k]*dV(i,j,k);
      }
    }    
  }

  /*--------------+
  |  x-direction  |
  +--------------*/
  CIPCSLx1(scheme.f    ,scheme.sigx);
  CIPCSLx2(scheme.sigy ,sxyz       );  /* sxy */
  CIPCSLx3(scheme.sigz ,sxyz       );  /* szx */
  CIPCSLx4(sxyz        ,clr        );  /* syz */

  /*---------------------+
  |  boundary condition  |
  +---------------------*/
  /* invert rho to clr */
  for(int i=ist; i<=ied; i++){
    for(int j=jst; j<=jed; j++){
      for(int k=kst; k<=ked; k++){
        clr[i][j][k]=clr[i][j][k]/dV(i,j,k);
      }
    }    
  }

  //clr.bnd_update();
  bdcond(clr); 
  clr.exchange_all();
  scheme.bdcond_f(clr);
  scheme.bdcond_i(clr);
  scheme.bdcond_j(clr);
  scheme.bdcond_k(clr);
  bdphiface(sxyz,Comp::i(),clr);
  bdphiface(sxyz,Comp::j(),clr);
  bdphiface(sxyz,Comp::k(),clr);
#ifdef IB
  ib_ext_scalar(clr);
  ib_bdcond(clr);
#endif
  scheme.f.exchange_all();
  scheme.sigx.exchange_all();
  scheme.sigy.exchange_all();
  scheme.sigz.exchange_all();
  sxyz.exchange_all();

#if 1
  /* convert clr to rho */
  for(int i=ist; i<=ied; i++){
    for(int j=jst; j<=jed; j++){
      for(int k=kst; k<=ked; k++){
        clr[i][j][k]=clr[i][j][k]*dV(i,j,k);
      }
    }    
  }

  /*--------------+
  |  y-direction  |
  +--------------*/
  CIPCSLy1(scheme.f    ,scheme.sigy);
  CIPCSLy2(scheme.sigx ,sxyz       ); //sxy
  CIPCSLy3(scheme.sigz ,sxyz       ); //syz
  CIPCSLy4(sxyz        ,clr        ); //szx
  //debug:bdphiface_check(sxyz,clr);

  /*---------------------+
  |  boundary condition  |
  +---------------------*/
  /* invert rho to clr */
  for(int i=ist; i<=ied; i++){
    for(int j=jst; j<=jed; j++){
      for(int k=kst; k<=ked; k++){
        clr[i][j][k]=clr[i][j][k]/dV(i,j,k);
      }
    }    
  }

  //clr.bnd_update();
  bdcond(clr); 
  clr.exchange_all();
  scheme.bdcond_f(clr);
  scheme.bdcond_i(clr);
  scheme.bdcond_j(clr);
  scheme.bdcond_k(clr);
  bdphiface(sxyz,Comp::i(),clr);
  bdphiface(sxyz,Comp::j(),clr);
  bdphiface(sxyz,Comp::k(),clr);
#ifdef IB
  ib_ext_scalar(clr);
  ib_bdcond(clr);
#endif
  scheme.f.exchange_all();
  scheme.sigx.exchange_all();
  scheme.sigy.exchange_all();
  scheme.sigz.exchange_all();
  sxyz.exchange_all();
#endif

  /* convert clr to rho */
  for(int i=ist; i<=ied; i++){
    for(int j=jst; j<=jed; j++){
      for(int k=kst; k<=ked; k++){
        clr[i][j][k]=clr[i][j][k]*dV(i,j,k);
      }
    }    
  }

  /*--------------+
  |  z-direction  |
  +--------------*/
  CIPCSLz1(scheme.f    ,scheme.sigz);
  CIPCSLz2(scheme.sigx ,sxyz       ); //szx
  CIPCSLz3(scheme.sigy ,sxyz       ); //syz
  CIPCSLz4(sxyz        ,clr        ); //sxy

  /*---------------------+
  |  boundary condition  |
  +---------------------*/
  /* invert rho to clr */
  for(int i=ist; i<=ied; i++){
    for(int j=jst; j<=jed; j++){
      for(int k=kst; k<=ked; k++){
        clr[i][j][k]=clr[i][j][k]/dV(i,j,k);
      }
    }    
  }

  //clr.bnd_update();
  bdcond(clr); 
  clr.exchange_all();
  scheme.bdcond_f(clr);
  scheme.bdcond_i(clr);
  scheme.bdcond_j(clr);
  scheme.bdcond_k(clr);
  bdphiface(sxyz,Comp::i(),clr);
  bdphiface(sxyz,Comp::j(),clr);
  bdphiface(sxyz,Comp::k(),clr);
#ifdef IB
  ib_ext_scalar(clr);
  ib_bdcond(clr);
#endif
  scheme.f.exchange_all();
  scheme.sigx.exchange_all();
  scheme.sigy.exchange_all();
  scheme.sigz.exchange_all();
  sxyz.exchange_all();

#ifdef IB
  ib_ext_scalar(clr);
#endif

#if 0
  if(time->current_step()==1){
    plot_f("f");
    plot_sigx("sigx");
    plot_sigy("sigy");
    plot_sigz("sigz");
    plot_sxyz("sxyzi",Comp::i());
    plot_sxyz("sxyzj",Comp::j());
    plot_sxyz("sxyzk",Comp::k());
    exit(0);
  }
#endif

  boil::cart.sum_real(&sum_outlet);
  boil::cart.sum_real(&sum_outletm);
  boil::timer.stop("cipcsl2 convection");

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_convection.cpp,v 1.12 2015/06/29 18:24:29 sato Exp $'/
+-----------------------------------------------------------------------------*/
