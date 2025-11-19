#include "nucleation.h"

/***************************************************************************//**
* set loop range
*******************************************************************************/
void Nucleation::set_range(std::vector<Site> & s) {

  real rng = rseed + eps;
  int ns=s.size()-1;
  real xs = s[ns].x()-rng;
  real xe = s[ns].x()+rng;
  real ys = s[ns].y()-rng;
  real ye = s[ns].y()+rng;
  real zs = s[ns].z()-rng;
  real ze = s[ns].z()+rng;
  //std::cout<<ns<<" "<<xs<<" "<<xe<<" "<<ys<<" "<<ye<<" "<<zs<<" "<<ze<<"\n";

  int ic = vf->aim(s[ns].x(), boil::pico);
  int jc = vf->ajm(s[ns].y(), boil::pico);
  int kc = vf->akm(s[ns].z(), boil::pico);

  int is = vf->aim(xs, boil::pico);
  int ie = vf->aip(xe, boil::pico);
  int js = vf->ajm(ys, boil::pico);
  int je = vf->ajp(ye, boil::pico);
  int ks = vf->akm(zs, boil::pico);
  int ke = vf->akp(ze, boil::pico);
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

  return;
}
