#include "vof.h"

/* liq+gas, both upwind */
real VOF::calc_diabatic_flux(const real jv, const real gliq, const real ggas,
                             real phiup,
                             const real nx, const real ny, const real nz) {
  if(fabs(gliq)<boil::pico) {
    return 0.0;
  } else {
    assert(gliq*ggas>=0.0);
    real cflrat = ggas/gliq; 

    real vn1 = -nx;
    real vn2 = -ny;
    real vn3 = -nz;
 
    real vm1 = fabs(vn1);
    real vm2 = fabs(vn2);
    real vm3 = fabs(vn3)+boil::pico;
    real qa = 1.0/(vm1+vm2+vm3);
    vm1 *= qa;
    vm2 *= qa;
    vm3 *= qa;

    real alphaliq = calc_alpha(phiup, vm1, vm2, vm3);
    real phigas = 1.-phiup;
    real alphagas = calc_alpha(phigas,vm1,vm2,vm3);

    real sgnj = (jv>0.0) - (jv<0.0);
    real dirsgn = (vn1*gliq>0.0) - (vn1*gliq<0.0);
#if 1
    return sgnj*iterate_flux(sgnj*jv,dirsgn,cflrat,
                             alphaliq,alphagas,vm1,vm2,vm3);
#endif
  }

  return 0.0;
}

/* liq upw, gas downwind */
real VOF::calc_diabatic_flux(real & jv, const real gliq, const real ggas,
                             const real dxrat, real phiup, real phidn,
                             const real nxup, const real nyup, const real nzup,
                             const real nxdn, const real nydn, const real nzdn) {

  if(fabs(gliq)<boil::pico) {
    return 0.0;
  } else {
    assert(gliq*ggas<=0.0);
    real cflrat = -ggas/gliq*dxrat; 

    real vn1up = -nxup;
    real vn2up = -nyup;
    real vn3up = -nzup;
 
    real vm1up = fabs(vn1up);
    real vm2up = fabs(vn2up);
    real vm3up = fabs(vn3up)+boil::pico;
    real qaup = 1.0/(vm1up+vm2up+vm3up);
    vm1up *= qaup;
    vm2up *= qaup;
    vm3up *= qaup;

    real alphaup = calc_alpha(phiup, vm1up, vm2up, vm3up);

    real vn1dn = -nxdn;
    real vn2dn = -nydn;
    real vn3dn = -nzdn;
 
    real vm1dn = fabs(vn1dn);
    real vm2dn = fabs(vn2dn);
    real vm3dn = fabs(vn3dn)+boil::pico;
    real qadn = 1.0/(vm1dn+vm2dn+vm3dn);
    vm1dn *= qadn;
    vm2dn *= qadn;
    vm3dn *= qadn;

    real phigas = 1.-phidn;
    real alphadn = calc_alpha(phigas, vm1dn, vm2dn, vm3dn);

    real dirsgnup = (vn1up*gliq>0.0) - (vn1up*gliq<0.0);
    real dirsgndn = (vn1dn*ggas>0.0) - (vn1dn*ggas<0.0);
 
    /* can gas cross the boundary in the geometrical sense? */
    /* gas flows in the opposite direction and sense of n is inverted,
     * dirsgn -> -dirsgn and minus in front */
    real testflux = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,boil::micro,-dirsgndn);
    /* it cannot cross */
    if(fabs(testflux)<boil::pico) {
      /* unphysical case */
      if(jv*gliq<=0.0) {
        return 0.0;
      } else {
        return jv;
      } 
    } else {
       /* is j aligned with liquid flux? */
       real sgng = (gliq>0.0) - (gliq<0.0); 
       if(jv*gliq>=0.0) {
         jv = fabs(jv);
       } else {
         jv = -fabs(jv);
       }
       real maxliq = calc_v_iter(alphaup,vm1up,vm2up,vm3up,
                                 flux_cfl,dirsgnup);
       real maxgasx = std::min(1.0,flux_cfl*cflrat);
       real maxgas = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,
                                  maxgasx,-dirsgndn);
       boil::aout<<"diabatic,case3 "<<maxliq<<" "<<maxgas<<" "<<maxliq+maxgas-jv<<boil::endl;
       /* do 0.0 and flux_cfl provide sufficient bounds? */
       if((maxliq+maxgas-jv)*jv>0.0) {
#if 1
         return sgng*iterate_flux(jv,dirsgnup,dirsgndn,cflrat,
                                  alphaup,vm1up,vm2up,vm3up,
                                  alphadn,vm1dn,vm2dn,vm3dn,
                                  0.0,flux_cfl);
#endif
       } else {
         real testx = 0.0;
         real testgasx;
         real testval = -jv;
         real testliq, testgas;
         real lastx, lastval;
         real loop(0.0);
         real upper_limit = real(maxiter);
         bool converged(false);
         while(!converged) {
           lastx = testx;
           lastval = testval;
           loop += 1.0;
           if(loop>upper_limit-0.5) {
             break;
           } 
           testx = loop*flux_cfl/upper_limit;
           testliq = calc_v_iter(alphaup,vm1up,vm2up,vm3up,
                                testx,dirsgnup);
           testgasx = std::min(1.0,testx*cflrat);
           testgas = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,
                                  testgasx,-dirsgndn);
           testval = testliq+testgas-jv;
           if       (testval*lastval<=0.0) {
             converged = true;
           }
         } /* while loop */
         if(!converged) {
           return 0.0;
         } else {
           boil::aout<<"diabatic,case3,iter "<<lastx<<" "<<testx<<" "<<lastval<<" "<<testval<<boil::endl;
#if 1
           return sgng*iterate_flux(jv,dirsgnup,dirsgndn,cflrat,
                                    alphaup,vm1up,vm2up,vm3up,
                                    alphadn,vm1dn,vm2dn,vm3dn,
                                    lastx,testx);
#endif

         }
       } /* sgn at 0.0 and flux_cfl is the same */
    } /* gas can cross boundary */
  } /* there is liquid flow */

  return 0.0;
}

/* underdevelopment */
#if 0
/* liq+gas upw and gas downwind */
real VOF::calc_diabatic_flux(const real jv, const real gliq, 
                             const real ggasup,const real ggasdn,
                             const real dxrat, real phiup, real phidn,
                             const real nxup, const real nyup, const real nzup,
                             const real nxdn, const real nydn, const real nzdn) {
  if(fabs(gliq)<boil::pico) {
    return 0.0;
  } else {
    assert(gliq*ggasup>=0.0);
    assert(gliq*ggasdn<=0.0);
    real cflratup = ggasup/gliq; 
    real cflratdn = ggasdn/gliq; 

    real vn1up = -nxup;
    real vn2up = -nyup;
    real vn3up = -nzup;
 
    real vm1up = fabs(vn1up);
    real vm2up = fabs(vn2up);
    real vm3up = fabs(vn3up)+boil::pico;
    real qaup = 1.0/(vm1up+vm2up+vm3up);
    vm1up *= qaup;
    vm2up *= qaup;
    vm3up *= qaup;

    real alphaupliq = calc_alpha(phiup, vm1up, vm2up, vm3up);
    real phiupgas = 1.-phiup;
    real alphaupgas = calc_alpha(phiupgas,vm1up,vm2up,vm3up);

    real vn1dn = -nxdn;
    real vn2dn = -nydn;
    real vn3dn = -nzdn;
 
    real vm1dn = fabs(vn1dn);
    real vm2dn = fabs(vn2dn);
    real vm3dn = fabs(vn3dn)+boil::pico;
    real qadn = 1.0/(vm1dn+vm2dn+vm3dn);
    vm1dn *= qadn;
    vm2dn *= qadn;
    vm3dn *= qadn;

    real phidngas = 1.-phidn;
    real alphadn = calc_alpha(phidngas, vm1dn, vm2dn, vm3dn);

    real sgnj = (jv>0.0) - (jv<0.0);
#if 0
    return sgnj*iterate_flux(sgnj*jv,cflratup,cflratdn,dxrat,
                             alphaupliq,alphaupgas,vm1up,vm2up,vm3up,
                                           alphadn,vm1dn,vm2dn,vm3dn);
#endif
  }

  return 0.0;
}
#endif
