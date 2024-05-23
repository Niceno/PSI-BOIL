#include "nucleation.h"
using namespace std;

/***************************************************************************//**
* set loop range
*******************************************************************************/
void Nucleation::set_range(std::vector<Site> & s) {

  /* define ns */
  int ns=s.size()-1;

  /* set_contain */
  /* center of bubble is inside/outside of decomposed domain */
  s[ns].set_contain_cb(vf->domain()->contains_xyz(s[ns].x(),s[ns].y(),s[ns].z()));
  /* nucleation site (assumed to be z=0) is inside/outside of decomposed domain */
  s[ns].set_contain_site(vf->domain()->contains_xyz(s[ns].x(),s[ns].y(),-boil::nano));

  /* set range for loop around a nucleation site */
  real rng = rseed + eps;
  real xs = s[ns].x()-rng;
  real xe = s[ns].x()+rng;
  real ys = s[ns].y()-rng;
  real ye = s[ns].y()+rng;
  real zs = s[ns].z()-rng;
  real ze = s[ns].z()+rng;
  //std::cout<<"proc= "<<boil::cart.iam()<<" ns= "<<ns<<" x "<<xs<<" "<<xe
  //<<" y "<<ys<<" "<<ye<<" z "<<zs<<" "<<ze<<"\n";

  /* set_contain_range: range is inside/outside of decomposed domain */
  /* Note: <= and >= are used for the detection */
  int sum_contains = vf->domain()->contains_xyz( xs, ys, zs)
                   + vf->domain()->contains_xyz( xs, ys, ze)
                   + vf->domain()->contains_xyz( xs, ye, zs)
                   + vf->domain()->contains_xyz( xs, ye, ze)
                   + vf->domain()->contains_xyz( xe, ys, zs)
                   + vf->domain()->contains_xyz( xe, ys, ze)
                   + vf->domain()->contains_xyz( xe, ye, zs)
                   + vf->domain()->contains_xyz( xe, ye, ze);
  //std::cout<<"proc= "<<boil::cart.iam()<<" sum_contains= "<<sum_contains<<"\n";
  s[ns].set_contain_range(false);
  if (sum_contains >=1) s[ns].set_contain_range(true);

  if (s[ns].contain_site()) {
    /* WARNING!!! ic, jc, kc are defined using minus (aim, ajm, akm) */
    //int ic = vf->aim(s[ns].x(), boil::femto);
    //int jc = vf->ajm(s[ns].y(), boil::femto);
    //int kc = vf->akm(s[ns].z(), boil::femto);
    int icw = vf->i(s[ns].x());
    int jcw = vf->j(s[ns].y());
    int kcw = vf->k(-boil::nano);
    //std::cout<<"set_range: proc= "<<boil::cart.iam()<<" ns= "<<ns<<" contain "<<s[ns].contain()
    //  <<" contain_range "<<s[ns].contain_range()<<" "
    //  <<ic<<" "<<jc<<" "<<kc<<" "<<s[ns].z()<<"\n";
    //boil::oout<<vf->xn(19)<<" "<<vf->xn(20)<<"\n";
    //exit(0);

    s[ns].set_ic_site(icw);
    s[ns].set_jc_site(jcw);
    s[ns].set_kc_site(kcw);
  }

  if (s[ns].contain_range()) {
    //int is = vf->aim(xs, boil::femto);
    //int ie = vf->aip(xe, boil::femto);
    //int js = vf->ajm(ys, boil::femto);
    //int je = vf->ajp(ye, boil::femto);
    //int ks = vf->akm(zs, boil::femto);
    //int ke = vf->akp(ze, boil::femto);
    int is = max(vf->aim(xs, boil::femto),vf->si());
    int ie = min(vf->aip(xe, boil::femto),vf->ei());
    int js = max(vf->ajm(ys, boil::femto),vf->sj());
    int je = min(vf->ajp(ye, boil::femto),vf->ej());
    int ks = max(vf->akm(zs, boil::femto),vf->sk());
    int ke = min(vf->akp(ze, boil::femto),vf->ek());
    //std::cout<<"set_range: proc= "<<boil::cart.iam()<<" "<<ns<<" "<<is<<" "<<ie<<" "
    //         <<js<<" "<<je<<" "<<ks<<" "<<ke<<"\n";

    s[ns].set_is(is);
    s[ns].set_ie(ie);
    s[ns].set_js(js);
    s[ns].set_je(je);
    s[ns].set_ks(ks);
    s[ns].set_ke(ke);
  }

  return;
}
