#include "vof.h"

/******************************************************************************/
void VOF::advance_x() {
  
  /* advance in the x-direction */

  Comp m = Comp::u();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei()+1; i++)
  //for(int j = sj(); j <=ej(); j++)
  //for(int k = sk(); k <=ek(); k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;
     
    real jv   = (*u)[m][i][j][k];
    real uval = vel_value(m,i,j,k);
    /* upwind i-index */
    int iup(i-1), idn(i);
    if(uval<0.0) {
      iup = i; 
      idn = i-1;
    }            
    real dxup = phi.dxc(iup);
    real dt = time->dt();

    /* fext is related to mflx as follows:
    * fext = -m'''/rhol = -m''*adens/rhol */
    real fextup = fext[iup][j][k];
    real fextdn = fext[idn][j][k];

    /* volumetric conservation factor, calculated in vof_advance
     * this factor represents the renormalisation of phi' after pc and must
     * be thus included in the real space flux after iteration */
    real volfactorup = stmp3[iup][j][k];
    real volfactordn = stmp3[idn][j][k];

    bool sourceup(true), sourcedn(true);

    if(fabs(fextup)<boil::pico) sourceup = false; 
    if(fabs(fextdn)<boil::pico) sourcedn = false; 

    real phiup = phi[iup][j][k];
    real phidn = phi[idn][j][k];

    if       (dxup==0.0||phiup<boil::pico) {
      if(jv*uval<0.0) {
        f = 0.0;
      } else {
        f = phiup * dSx(i,j,k) * jv * dt;
      }
    } else if(phiup>(1.0-boil::pico)) {
      if       (jv*uval<0.0) {
        f = 0.0;
      } else if(phidn>(1.0-boil::pico)) {
        f = phiup * dSx(i,j,k) * jv * dt;
#if 0
        if(fabs(f)<boil::micro) {
          f = 0.0; /* remove spurious fluxes */
        }
#endif
      } else {
        f = phiup * dSx(i,j,k) * jv * dt;
      }    
    //} else if(!sourceup&!sourcedn) {
    } else if(true) {
      /* calculate g: CFL upwind */
      real g = uval*dt/dxup;
      f = dV(iup,j,k)*calc_flux(g,phiup,nx[iup][j][k],
                                        ny[iup][j][k],
                                        nz[iup][j][k]);
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
        real adensup = adens[iup][j][k]+boil::pico;
        real mmx = mx[iup][j][k]; /* mx towards the liquid in real space */
        gasvelup = uval + (fextup*denscoef/adensup)*(-mmx);

        /* does gas flow in the same direction as liquid? */
        if(gasvelup*uval>=0.0)
          gasflowup++;
      }
      if(sourcedn) {
        gasflowdn++;
        real adensdn = adens[idn][j][k]+boil::pico;
        real mmx = mx[idn][j][k];
        gasveldn = uval + (fextdn*denscoef/adensdn)*(-mmx);
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

      //boil::aout<<i<<" "<<j<<" "<<k<<" | "<<iup<<" "<<sourceup<<" "<<sourcedn<<" | "<<jv*dt/dxup<<" "<<uval<<" "<<gasvelup<<" "<<gasveldn<<" | "<<gasflowup<<" "<<gasflowdn<<boil::endl;

      /* case 1: downwind follows flow, upwind opposes: only liq crosses bnd */
      if       (gasflowup==0&&gasflowdn==1) {
        /* unphysical case */
        if(jv*uval<0.0) {
          f = 0.0;
        } else {
          f = dSx(i,j,k) * jv * dt;
        }
        boil::aout<<"case1x "<<iup<<" "<<j<<" "<<k<<" | "<<idn<<" "<<phiup<<" "<<stmp[iup][j][k]/dV(iup,j,k)<<" "<<f/dV(iup,j,k)<<" | "<<(stmp[idn][j][k]+fabs(f))/dV(idn,j,k)<<" "<<(stmp[iup][j][k]-fabs(f))/dV(iup,j,k)<<boil::endl;
      /* case 2: both downwind and upwind follow the flow */
      } else if(gasflowup==1||gasflowdn==1) { 
        /* calculate gliq and ggas: CFL of liquid and gas upwind */
        real gliq = uval*dt/dxup;
        real ggas;
        if(gasflowup==1) {
          ggas = gasvelup*dt/dxup;
        } else {
          ggas = gasveldn*dt/dxup;
        }
        /* dimensionless vol. flow = jv*dt*dS/dV*volfactor */
        real gj = jv*dt/dxup/volfactorup; 
#if 1
        f = dV(iup,j,k)*volfactorup
            *calc_diabatic_flux(gj,gliq,ggas,phiup,
                                nx[iup][j][k],ny[iup][j][k],nz[iup][j][k]);
        boil::aout<<"case2x "<<iup<<" "<<j<<" "<<k<<" | "<<idn<<" "<<phiup<<" "<<stmp[iup][j][k]/dV(iup,j,k)<<" "<<f/dV(iup,j,k)<<" | "<<(stmp[idn][j][k]+fabs(f))/dV(idn,j,k)<<" "<<(stmp[iup][j][k]-fabs(f))/dV(iup,j,k)<<" "<<gliq<<" "<<ggas<<" | "<<gj<<" "<<uval<<" "<<(*u)[m][i][j][k]<<" "<<uliq[m][i][j][k]<<boil::endl;
#endif

      /* case 3: both downwind and upwind oppose the flow */ 
      } else if(gasflowdn==0||gasflowup==0) {
        real dxdn = phi.dxc(idn);

        /* degenerate case */
        if(dxdn==0.0) {
          real jgas;
          if(gasflowdn==0) {
            jgas = (1.0-phi[idn][j][k])*gasveldn;
          } else {
            jgas = (1.0-phi[idn][j][k])*gasvelup;
          }
          jv -= jgas;

          /* unphysical case */
          if(jv*uval<0.0) {
            f = 0.0;
          } else {
            f = dSx(i,j,k)*jv*dt;
          } 
          boil::aout<<"case3x "<<iup<<" "<<j<<" "<<k<<" | "<<idn<<" "<<phiup<<" "<<stmp[iup][j][k]/dV(iup,j,k)<<" "<<f/dV(iup,j,k)<<" | "<<(stmp[idn][j][k]+fabs(f))/dV(idn,j,k)<<" "<<(stmp[iup][j][k]-fabs(f))/dV(iup,j,k)<<boil::endl;
        } else {
          /* calculate gliq and ggas: CFL of liquid upw and gas dnw */
          real gliq = uval*dt/dxup;
          real ggas;
          if(gasflowdn==0) {
            ggas = gasveldn*dt/dxdn;
          } else {
            ggas = gasvelup*dt/dxdn;
          }
          real gj = jv*dt/dxup/volfactorup;

          /* calculate ratio of dx, including the effect of vol expansion */
          real dxrat = dxdn/dxup*volfactordn/volfactorup;
#if 1
          f = volfactorup*dV(iup,j,k)
              *calc_diabatic_flux(gj,gliq,ggas,dxrat,phiup,phidn,
                                  nx[iup][j][k],ny[iup][j][k],nz[iup][j][k],
                                  nx[idn][j][k],ny[idn][j][k],nz[idn][j][k]);
          boil::aout<<"case3x "<<iup<<" "<<j<<" "<<k<<" | "<<idn<<" "<<phiup<<" "<<stmp[iup][j][k]/dV(iup,j,k)<<" "<<f/dV(iup,j,k)<<" | "<<(stmp[idn][j][k]+fabs(f))/dV(idn,j,k)<<" "<<(stmp[iup][j][k]-fabs(f))/dV(iup,j,k)<<boil::endl;
#endif
        }
      /* case 4: downwind opposes flow, upwind follows: all three cross bnd */
      } else {
        /* underdevelopment */
        boil::oout<<"VOF_advance_x: underdevelopment! Exiting."<<boil::endl;
        exit(0);
#if 0
        real dxdn = phi.dxc(idn);

        /* degenerate case */
        if(dxdn==0.0) {
          real jgas = (1.0-phi[idn][j][k])*gasveldn;
          jv -= jgas;

          /* unphysical case */
          if(jv*uval<0.0) {
            f = 0.0;
          } else {
            real gliq = uval*dt/dxup;
            real ggas = gasvelup*dt/dxup;
            real gj = jv*dt/dxup/volfactorup;
  #if 1
            f = volfactorup*dV(iup,j,k)
                *calc_diabatic_flux(gj,gliq,ggas,phiup,
                                    nx[iup][j][k],ny[iup][j][k],nz[iup][j][k]);
  #endif
          }
        } else {
          /* calculate gliq and ggas: CFL of liquid and gas upwind */
          real gliq = uval*dt/dxup;
          real ggasup = gasvelup*dt/dxup;
          real ggasdn = gasveldn*dt/dxdn;
          real gj = jv*dt/dxup;

          /* calculate ratio of dx, including the effect of vol expansion */
          real dxrat = dxdn/dxup*volfactordn/volfactorup;

          real phidn = phi[idn][j][k];
  #if 1
          f = volfactorup*dV(iup,j,k)
              *calc_diabatic_flux(gj,gliq,ggasup,ggasdn,dxrat,phiup,phidn,
                                  nx[iup][j][k],ny[iup][j][k],nz[iup][j][k],
                                  nx[idn][j][k],ny[idn][j][k],nz[idn][j][k]);
  #endif
        }
#endif
      } /* cases */
#endif
    } /* adiabatic-diabatic */

#if 0
    /* limit flux */
    real sgnf = (f>0.0)-(f<0.0);
    f = sgnf*std::min(fabs(stmp[iup][j][k]),fabs(f));
#endif

#if 0
    /* update stmp */
    stmp[i-1][j][k] = stmp[i-1][j][k] - f;
    stmp[i  ][j][k] = stmp[i  ][j][k] + f;
#else
    fluxmax[m][i][j][k] = f;
    sosflux[m][i][j][k] = f/3.0;
#endif

#if 0
    if((i==100||i==101)&&j==3&&k==100) {
      std::cout<<"advance_x:"<<f<<" "<<i<<" "<<phi[i-1][j][k]<<"\n";
      //if (phi[i-1][j][k]>0.1) exit(0);
    }
#endif

  }

  return;

}

