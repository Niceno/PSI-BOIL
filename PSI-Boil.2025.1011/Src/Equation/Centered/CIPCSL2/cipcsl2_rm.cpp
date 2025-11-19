#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::rm(char * nm, const int it) {
/***************************************************************************//**
*  \brief Save variables of CIPCSL
*******************************************************************************/
  std::string names;

  /* phi */
  names = nm + std::string("-phi");
  char *cphi;
  cphi=&names[0];
  phi.rm(cphi,it);

  /* clr */
  names = nm + std::string("-clr");
  char *ccell;
  ccell=&names[0];
  clr.rm(ccell,it);

  /* node */
  names = nm + std::string("-f");
  char *cnode;
  cnode=&names[0];
  scheme.f.rm(cnode,it);

  /* edge */
  names = nm + std::string("-sigx");
  char *csigx;
  csigx=&names[0];
  scheme.sigx.rm(csigx,it);

  names = nm + std::string("-sigy");
  char *csigy;
  csigy=&names[0];
  scheme.sigy.rm(csigy,it);

  names = nm + std::string("-sigz");
  char *csigz;
  csigz=&names[0];
  scheme.sigz.rm(csigz,it);

  /* face */
  names = nm + std::string("-face");
  char *cface;
  cface=&names[0];
  sxyz.rm(cface,it);
}
