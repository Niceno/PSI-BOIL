#include "nucleation.h"

/***************************************************************************//**
* set loop range
*******************************************************************************/
void Nucleation::set_range(std::vector<Site> & s) {

  /* define ns */
  int ns=s.size()-1;

  /* set_contain */
  /* nucleation site (assumed to be z=0) is inside/outside of decomposed domain */
  s[ns].set_contain_site(clr->domain()->contains_xyz(s[ns].x(),s[ns].y(),-boil::nano));

  /* set range for loop around a nucleation site */
  real eps = 1.0 * dxmin;
  real rng = rseed + eps;
  real xs = s[ns].x()-rng;
  real xe = s[ns].x()+rng;
  real ys = s[ns].y()-rng;
  real ye = s[ns].y()+rng;
  real zs = s[ns].z()-rng;
  real ze = s[ns].z()+rng;
#if 0
  if (ns==349) {
    boil::oout<<"set_range:ns="<<ns<<" "<<xs<<" "<<xe<<" "<<ys<<" "<<ye<<" "<<zs<<" "<<ze<<"\n";
  }
#endif

  /* set_contain_range: range is inside/outside of decomposed domain */
  /* Note: <= and >= are used for the detection */
  int sum_contains = clr->domain()->contains_xyz( xs, ys, zs)
                   + clr->domain()->contains_xyz( xs, ys, ze)
                   + clr->domain()->contains_xyz( xs, ye, zs)
                   + clr->domain()->contains_xyz( xs, ye, ze)
                   + clr->domain()->contains_xyz( xe, ys, zs)
                   + clr->domain()->contains_xyz( xe, ys, ze)
                   + clr->domain()->contains_xyz( xe, ye, zs)
                   + clr->domain()->contains_xyz( xe, ye, ze);
  //std::cout<<"set_range:ns= "<<ns<<" sum_contains= "<<sum_contains
  //         <<" x "<<s[ns].x()<<" y "<<s[ns].y()<<"\n";

  s[ns].set_contain_range(false);
  if (sum_contains >=1) s[ns].set_contain_range(true);

  /* set (i,j,k) of nucleation site in wall */
  if (s[ns].contain_site()) {
    int icw = clr->i(s[ns].x());
    int jcw = clr->j(s[ns].y());
    int kcw = clr->k(-boil::nano); // in solid!
    s[ns].set_ic_site(icw);
    s[ns].set_jc_site(jcw);
    s[ns].set_kc_site(kcw);
  }

  if (s[ns].contain_range()) {
    int is = clr->aim(xs, boil::femto);
    int ie = clr->aip(xe, boil::femto);
    int js = clr->ajm(ys, boil::femto);
    int je = clr->ajp(ye, boil::femto);
    int ks = clr->akm(zs, boil::femto);
    int ke = clr->akp(ze, boil::femto);
    //std::cout<<ns<<" "<<ic<<" "<<jc<<" "<<kc<<"\n";
    //std::cout<<ns<<" "<<is<<" "<<ie<<" "<<js<<" "<<je<<" "<<ks<<" "<<ke<<"\n";
    is = std::max(is,clr->si());
    ie = std::min(ie,clr->ei());
    js = std::max(js,clr->sj());
    je = std::min(je,clr->ej());
    ks = std::max(ks,clr->sk());
    ke = std::min(ke,clr->ek());

    s[ns].set_is(is);
    s[ns].set_ie(ie);
    s[ns].set_js(js);
    s[ns].set_je(je);
    s[ns].set_ks(ks);
    s[ns].set_ke(ke);
    //std::cout<<"set_range:ns= "<<ns<<" sum_contains= "<<sum_contains
    //       <<" x "<<s[ns].x()<<" y "<<s[ns].y()
    //       <<" isie "<<s[ns].is()<<" "<<s[ns].ie()
    //       <<" jsje "<<s[ns].js()<<" "<<s[ns].je()<<"\n";

    /* kadj: k in adjacent cell of wall (fluid side) */
    int kadj = 0;
    if (s[ns].contain_range() ) {
      for(int k=s[ns].ks(); k<=s[ns].ke(); k++) {
        int i = s[ns].is();
        int j = j=s[ns].js();
        if(clr->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
          if(clr->domain()->ibody().off(i,j,k-1)) {
            kadj=k;
          }
        }
      }
    }
    s[ns].set_kadj(kadj);
  }

  return;
}
