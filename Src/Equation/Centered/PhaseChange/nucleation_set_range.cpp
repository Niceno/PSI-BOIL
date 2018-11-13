#include "nucleation.h"

/***************************************************************************//**
* set loop range
*******************************************************************************/
void Nucleation::set_range(std::vector<Site> & s) {

  real eps = rseed + dxmin * 1.5;
  int ns=s.size()-1;
  real xs = s[ns].x()-eps;
  real xe = s[ns].x()+eps;
  real ys = s[ns].y()-eps;
  real ye = s[ns].y()+eps;
  real zs = s[ns].z()-eps;
  real ze = s[ns].z()+eps;
  //std::cout<<ns<<" "<<xs<<" "<<xe<<" "<<ys<<" "<<ye<<" "<<zs<<" "<<ze<<"\n";

  int ic = clr->aim(s[ns].x(), boil::pico);
  int jc = clr->ajm(s[ns].y(), boil::pico);
  int kc = clr->akm(s[ns].z(), boil::pico);

  int is = clr->aim(xs, boil::pico);
  int ie = clr->aip(xe, boil::pico);
  int js = clr->ajm(ys, boil::pico);
  int je = clr->ajp(ye, boil::pico);
  int ks = clr->akm(zs, boil::pico);
  int ke = clr->akp(ze, boil::pico);
  //std::cout<<ns<<" "<<ic<<" "<<jc<<" "<<kc<<"\n";
  //std::cout<<ns<<" "<<is<<" "<<ie<<" "<<js<<" "<<je<<" "<<ks<<" "<<ke<<"\n";

  s[ns].set_ic(ic);
  s[ns].set_jc(jc);
  s[ns].set_kc(kc);

  s[ns].set_is(is);
  s[ns].set_ie(ie);
  s[ns].set_js(js);
  s[ns].set_je(je);
  s[ns].set_ks(ks);
  s[ns].set_ke(ke);
}

