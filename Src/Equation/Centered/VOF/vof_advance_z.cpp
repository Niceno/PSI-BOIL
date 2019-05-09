#include "vof.h"

/******************************************************************************/
void VOF::advance_z() {

  /* advance in the z-direction */

  Comp m = Comp::w();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei(); i++)
  //for(int j = sj(); j <=ej(); j++)
  //for(int k = sk(); k <=ek()+1; k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;
     
    real jv   = (*u)[m][i][j][k];
    real uval = vel_value(m,i,j,k);
    /* upwind k-index */
    int kup(k-1), kdn(k);
    if(uval<0.0) {
      kup = k; 
      kdn = k-1;
    }            
    real dzup = phi.dzc(kup);
    real dt = time->dt();

    /* fext is related to mflx as follows:
    * fext = -m'''/rhol = -m''*adens/rhol */
    real fextup = fext[i][j][kup];
    real fextdn = fext[i][j][kdn];

    /* volumetric conservation factor, calculated in vof_advance
     * this factor represents the renormalisation of phi' after pc and must
     * be thus included in the real space flux after iteration */
    real volfactorup = stmp3[i][j][kup];
    real volfactordn = stmp3[i][j][kdn];

    bool sourceup(true), sourcedn(true);

    if(fabs(fextup)<boil::pico) sourceup = false; 
    if(fabs(fextdn)<boil::pico) sourcedn = false; 

    real phiup = phi[i][j][kup];
    real phidn = phi[i][j][kdn];

    if       (dzup==0.0||phiup<boil::pico) {
      if(jv*uval<0.0) {
        f = 0.0;
      } else {
        f = phiup * dSz(i,j,k) * jv * dt;
      }
    } else if(phiup>(1.0-boil::pico)) {
      if       (jv*uval<0.0) {
        f = 0.0;
      } else if(phidn>(1.0-boil::pico)) {
        f = phiup * dSz(i,j,k) * jv * dt;
#if 0
        if(fabs(f)<boil::micro) {
          f = 0.0; /* remove spurious fluxes */
        }
#endif
      } else {
        f = phiup * dSz(i,j,k) * jv * dt;
      }    
    //} else if(!sourceup&!sourcedn) {
    } else if(true) {
      /* calculate g: CFL upwind */
      real g = uval*dt/dzup;
      f = dV(i,j,kup)*calc_flux(g,phiup,nz[i][j][kup],
                                        ny[i][j][kup],
                                        nx[i][j][kup]);
      //f *= volfactorup;
    } else {
      exit(0);
#if 0
      /* gas velocities */
      real gasvelup(-1000.0),gasveldn(-1000.0);
      int gasflowup(-1),gasflowdn(-1);
      real denscoef = 1.0-rhol/rhov;

      if(sourceup) {
        gasflowup++;
        real adensup = adens[i][j][kup]+boil::pico;
        real mmz = mz[i][j][kup]; /* mz towards the liquid in real space */
        gasvelup = uval + (fextup*denscoef/adensup)*(-mmz);

        /* does gas flow in the same direction as liquid? */
        if(gasvelup*uval>=0.0)
          gasflowup++;
      }
      if(sourcedn) {
        gasflowdn++;
        real adensdn = adens[i][j][kdn]+boil::pico;
        real mmz = mz[i][j][kdn];
        gasveldn = uval + (fextdn*denscoef/adensdn)*(-mmz);
        if(gasveldn*uval>=0.0)
          gasflowdn++;
      }

      /* case 4: downwind opposes flow, upwind follows: all three cross bnd 
       * -> underdeveloped: treated as either case 2 or case 3 */
      if(gasflowup==1&&gasflowdn==0) {
        real gasveldiff = gasvelup+gasveldn;
        if(gasveldiff*uval<0.0) {
          gasflowup = 0;
          gasveldn = gasveldiff;
          gasvelup = gasveldiff;
        } else {
          gasflowdn = 1; 
          gasvelup = gasveldiff;
          gasveldn = gasveldiff;
        }
      }

      //boil::aout<<i<<" "<<j<<" "<<k<<" | "<<kup<<" "<<sourceup<<" "<<sourcedn<<" | "<<jv*dt/dzup<<" "<<uval<<" "<<gasvelup<<" "<<gasveldn<<" | "<<gasflowup<<" "<<gasflowdn<<boil::endl;

      /* case 1: downwind follows flow, upwind opposes: only liq crosses bnd */
      if       (gasflowup==0&&gasflowdn==1) {
        /* unphysical case */
        if(jv*uval<0.0) {
          f = 0.0;
        } else {
          f = dSz(i,j,k) * jv * dt;
        }
        boil::aout<<"case1z "<<i<<" "<<j<<" "<<kup<<" | "<<kdn<<" "<<phiup<<" "<<stmp[i][j][kup]/dV(i,j,kup)<<" "<<f/dV(i,j,kup)<<" | "<<(stmp[i][j][kdn]+fabs(f))/dV(i,j,kdn)<<" "<<(stmp[i][j][kup]-fabs(f))/dV(i,j,kup)<<boil::endl;
      /* case 2: both downwind and upwind follow the flow */
      } else if(gasflowup==1||gasflowdn==1) { 
        /* calculate gliq and ggas: CFL of liquid and gas upwind */
        real gliq = uval*dt/dzup;
        real ggas;
        if(gasflowup==1) {
          ggas = gasvelup*dt/dzup;
        } else {
          ggas = gasveldn*dt/dzup;
        }
        /* dimensionless vol. flow = jv*dt*dS/dV*volfactor */
        real gj = jv*dt/dzup/volfactorup; 
#if 1
        f = dV(i,j,kup)*volfactorup
            *calc_diabatic_flux(gj,gliq,ggas,phiup,
                                nz[i][j][kup],ny[i][j][kup],nx[i][j][kup]);
        boil::aout<<"case2z "<<i<<" "<<j<<" "<<kup<<" | "<<kdn<<" "<<phiup<<" "<<stmp[i][j][kup]/dV(i,j,kup)<<" "<<f/dV(i,j,kup)<<" | "<<(stmp[i][j][kdn]+fabs(f))/dV(i,j,kdn)<<" "<<(stmp[i][j][kup]-fabs(f))/dV(i,j,kup)<<boil::endl;
#endif

      /* case 3: both downwind and upwind oppose the flow */ 
      } else if(gasflowdn==0||gasflowup==0) {
        real dzdn = phi.dzc(kdn);

        /* degenerate case */
        if(dzdn==0.0) {
          real jgas;
          if(gasflowdn==0) {
            jgas = (1.0-phi[i][j][kdn])*gasveldn;
          } else {
            jgas = (1.0-phi[i][j][kdn])*gasvelup;
          }
          jv -= jgas;

          /* unphysical case */
          if(jv*uval<0.0) {
            f = 0.0;
          } else {
            f = dSz(i,j,k)*jv*dt;
          } 
          boil::aout<<"case3z "<<i<<" "<<j<<" "<<kup<<" | "<<kdn<<" "<<phiup<<" "<<stmp[i][j][kup]/dV(i,j,kup)<<" "<<f/dV(i,j,kup)<<" | "<<(stmp[i][j][kdn]+fabs(f))/dV(i,j,kdn)<<" "<<(stmp[i][j][kup]-fabs(f))/dV(i,j,kup)<<boil::endl;
        } else {
          /* calculate gliq and ggas: CFL of liquid upw and gas dnw */
          real gliq = uval*dt/dzup;
          real ggas;
          if(gasflowdn==0) {
            ggas = gasveldn*dt/dzdn;
          } else {
            ggas = gasvelup*dt/dzdn;
          }
          real gj = jv*dt/dzup/volfactorup;

          /* calculate ratio of dz, including the effect of vol expansion */
          real dzrat = dzdn/dzup*volfactordn/volfactorup;
#if 1
          f = volfactorup*dV(i,j,kup)
              *calc_diabatic_flux(gj,gliq,ggas,dzrat,phiup,phidn,
                                  nz[i][j][kup],ny[i][j][kup],nx[i][j][kup],
                                  nz[i][j][kdn],ny[i][j][kdn],nx[i][j][kdn]);
          boil::aout<<"case3z "<<i<<" "<<j<<" "<<kup<<" | "<<kdn<<" "<<phiup<<" "<<stmp[i][j][kup]/dV(i,j,kup)<<" "<<f/dV(i,j,kup)<<" | "<<(stmp[i][j][kdn]+fabs(f))/dV(i,j,kdn)<<" "<<(stmp[i][j][kup]-fabs(f))/dV(i,j,kup)<<boil::endl;
#endif
        }
      /* case 4: downwind opposes flow, upwind follows: all three cross bnd */
      } else {
        /* underdevelopment */
        boil::oout<<"VOF_advance_z: underdevelopment! Exiting."<<boil::endl;
        exit(0);
#if 0

#endif
      } /* cases */
#endif
    } /* adiabatic-diabatic */

#if 0
    /* limit flux */
    real sgnf = (f>0.0)-(f<0.0);
    f = sgnf*std::min(fabs(stmp[i][j][kup]),fabs(f));
#endif

#if 0
    /* update stmp */
    stmp[i][j][k-1] = stmp[i][j][k-1] - f;
    stmp[i][j][k  ] = stmp[i][j][k  ] + f;
#else
    fluxmax[m][i][j][k] = f;
    sosflux[m][i][j][k] = f/3.0;
#endif

#if 0
    if((k==100||k==101)&&j==3&&i==100) {
      std::cout<<"advance_z:"<<f<<" "<<k<<" "<<phi[i][j][k-1]<<"\n";
    }
#endif

  }

  return;

}

