#include "vofaxisym.h"

real VOFaxisym::calc_flux_axisymmetric(const real gg,
                                       const int i, const int j, const int k,
                                       const Comp & mcomp) {
  /* no advection */
  if(gg==0.0) {
    return 0.0;
  }

  /* cell dimensions and cut values */
#if 1
  real g, ratio; /* cut-cell-fraction in linear and volumetric terms, resp */
  real absg; /* and in abs value */
  real eta0; /* one of the edges of the cut-cell */
  real eta00; /* left edge of the cut-cell, must be identified and rescaled */
  real eta1 = phi.xn(i)/phi.dxc(i); /* left */
  real eta2 = phi.xn(i+1)/phi.dxc(i); /* right */
  if       (mcomp == Comp::w()) {
    g = gg;
    absg = fabs(g);
    ratio = g;
    eta0 = eta1;
    eta00 = eta0;
  } else if(mcomp == Comp::u()) {
    real etaf;
    /* in the Cartesian limit etaf->inf, g & ratio become gg */
    if(gg>0.0) {
      etaf = eta2;
#if 0
      g =  etaf - sqrt(etaf*etaf-2.*etaf*gg);
      //g =  etaf - sqrt(etaf*etaf-(eta1+eta2)*gg);
#else
      g = etaf*(1.-exp(-gg/etaf)); /* derived from cst velocity */
#endif
      absg = fabs(g);
      eta0 = etaf-g;
      if(!(absg>0.0)) {
        return 0.0;
      }
      eta00 = eta0/absg;
    } else {
      etaf = eta1;
#if 0
      g =  etaf - sqrt(etaf*etaf-2.*etaf*gg);
      //g =  etaf - sqrt(etaf*etaf-(eta1+eta2)*gg);
#else
      g = etaf*(1.-exp(-gg/etaf));
#endif
      absg = fabs(g);
      if(!(absg>0.0)) {
        return 0.0;
      }
      eta0 = etaf-g;
      eta00 = etaf/absg;
    }
    ratio = (etaf*etaf-eta0*eta0)/(eta2*eta2-eta1*eta1);
  } else {
    boil::aout<<"Error! Wrong component specified for flux calculations! "
              <<"Exiting."<<boil::endl;
    exit(0);
  }
#endif

  /* n points to the liquid */
  real vnx = -nx[i][j][k];
  real vnz = -nz[i][j][k];

  /* L1-normalized normal vector */
  real vmx = fabs(vnx);
  real vmz = fabs(vnz);
  real qa = 1.0/(vmx+vmz+boil::pico);
  vmx *= qa;
  vmz *= qa;
  vnx *= qa;
  vnz *= qa;
  
  real * vm1 = &vmx;
  real * vm2 = &vmz;
  real * vn1 = &vnx;
  real * vn2 = &vnz;
  if(mcomp == Comp::w()) {
    vm1 = &vmz;
    vm2 = &vmx;
    vn1 = &vnz;
    vn2 = &vnx;
  }
  
  /* normalized alpha value in the cell 
   * even if nalpha is not realistic, the checks at the beginning 
   * of calc_v_axisymmetric should handle it */
  real alpha = qa*nalpha[i][j][k];
      
  /* this number serves two purposes:
     - alpha offset due to origin shift
     - 1-ra is the new L1 norm of the vm vector */
  real ra = (*vm1) * (1.0 - absg);
  qa = 1.0/(1.0-ra);

  /* if directions are aligned, alpha offset has to be introduced
     if not, the mirroring of the coordinate system nullifies the offset  */
  if (g*(*vn1) > 0) alpha = alpha - ra;

  /* normal vector in direction of transport is rescaled due to
     augmentation of cell dimensions */
  *vm1 = (*vm1) * absg;
  *vn1 = (*vn1) * absg;

  /* calculate f: reduced flux (for unit volume)
     - the qa factor renormalizes the line equation
     - the ratio factor then scales the transported volume back to the org wedge 
     - it also imposes the correct sign */
  real Kdummy;
  real f = calc_v_axisymmetric(vnx*qa,alpha*qa,eta00,Kdummy) * ratio;
  //boil::oout<<i<<" "<<k<<" | "<<vnx*qa<<" "<<(vmx+vmz)*qa<<" "<<alpha*qa<<" | "<<eta00<<" "<<ratio<<" "<<f<<boil::endl;

  return f;
}

