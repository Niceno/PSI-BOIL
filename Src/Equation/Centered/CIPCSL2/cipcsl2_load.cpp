#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::load(const char * nm, const int it) {
/***************************************************************************//**
*  \brief Load variables of CIPCSL2
*******************************************************************************/
  std::string names;

  /* phi */
  names = nm + std::string("-phi");
  char *cphi;
  cphi=&names[0];
  phi.load(cphi,it);

  /* clr */
  names = nm + std::string("-clr");
  char *ccell;
  ccell=&names[0];
  clr.load(ccell,it);

  /* node */
  names = nm + std::string("-f");
  char *cnode;
  cnode=&names[0];
  scheme.f.load(cnode,it);

  /* edge */
  names = nm + std::string("-sigx");
  char *csigx;
  csigx=&names[0];
  scheme.sigx.load(csigx,it);

  names = nm + std::string("-sigy");
  char *csigy;
  csigy=&names[0];
  scheme.sigy.load(csigy,it);

  names = nm + std::string("-sigz");
  char *csigz;
  csigz=&names[0];
  scheme.sigz.load(csigz,it);

  /* face */
  names = nm + std::string("-face");
  char *cface;
  cface=&names[0];
  sxyz.load(cface,it);

#if 0
  plot_f("f.dat");
  plot_sigx("sigx.dat");
  plot_sigy("sigy.dat");
  plot_sigz("sigz.dat");
  plot_sxyz("sxyzi.dat",Comp::i());
  plot_sxyz("sxyzj.dat",Comp::j());
  plot_sxyz("sxyzk.dat",Comp::k());
#endif

}
