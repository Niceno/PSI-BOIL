#include "vof.h"

/******************************************************************************/
void VOF::advance_y() {

  /* advance in the y-direction */

  Comp m = Comp::v();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei(); i++)
  //for(int j = sj(); j <=ej()+1; j++)
  //for(int k = sk(); k <=ek(); k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;
     
    real jv   = (*u)[m][i][j][k];
    real uval = vel_value(m,i,j,k);
    /* upwind j-index */
    int jup(j-1), jdn(j);
    if(uval<0.0) {
      jup = j; 
      jdn = j-1;
    }            
    real dyup = phi.dyc(jup);
    real dt = time->dt();

    /* fext is related to mflx as follows:
    * fext = -m'''/rhol = -m''*adens/rhol */
    real fextup = fext[i][jup][k];
    real fextdn = fext[i][jdn][k];

    /* volumetric conservation factor, calculated in vof_advance
     * this factor represents the renormalisation of phi' after pc and must
     * be thus included in the real space flux after iteration */
    real volfactorup = stmp3[i][jup][k];
    real volfactordn = stmp3[i][jdn][k];

    bool sourceup(true), sourcedn(true);

    if(fabs(fextup)<boil::pico) sourceup = false; 
    if(fabs(fextdn)<boil::pico) sourcedn = false; 

    real phiup = phi[i][jup][k];
    real phidn = phi[i][jdn][k];

    if       (dyup==0.0||phiup<boil::pico) {
      if(jv*uval<0.0) {
        f = 0.0;
      } else {
        f = phiup * dSy(i,j,k) * jv * dt;
      }
    } else if(phiup>(1.0-boil::pico)) {
      if       (jv*uval<0.0) {
        f = 0.0;
      } else if(phidn>(1.0-boil::pico)) {
        f = phiup * dSy(i,j,k) * jv * dt;
#if 0
        if(fabs(f)<boil::micro) {
          f = 0.0; /* remove spurious fluxes */
        }
#endif
      } else {
        f = phiup * dSy(i,j,k) * jv * dt;
      }    
    //} else if(!sourceup&!sourcedn) {
    } else if(true) {
      /* calculate g: CFL upwind */
      real g = uval*dt/dyup;
      f = dV(i,jup,k)*calc_flux(g,phiup,ny[i][jup][k],
                                        nx[i][jup][k],
                                        nz[i][jup][k]);
    } else {
      /* gas velocities */
      real gasvelup(-1000.0),gasveldn(-1000.0);
      int gasflowup(-1),gasflowdn(-1);
      real denscoef = 1.0-rhol/rhov;

      if(sourceup) {
        gasflowup++;
        real adensup = adens[i][jup][k]+boil::pico;
        real mmy = my[i][jup][k]; /* my towards the liquid in real space */
        gasvelup = uval + (fextup*denscoef/adensup)*(-mmy);

        /* does gas flow in the same direction as liquid? */
        if(gasvelup*uval>=0.0)
          gasflowup++;
      }
      if(sourcedn) {
        gasflowdn++;
        real adensdn = adens[i][jdn][k]+boil::pico;
        real mmy = my[i][jdn][k];
        gasveldn = uval + (fextdn*denscoef/adensdn)*(-mmy);
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

      //boil::aout<<i<<" "<<j<<" "<<k<<" | "<<jup<<" "<<sourceup<<" "<<sourcedn<<" | "<<jv*dt/dyup<<" "<<uval<<" "<<gasvelup<<" "<<gasveldn<<" | "<<gasflowup<<" "<<gasflowdn<<boil::endl;

      /* case 1: downwind follows flow, upwind opposes: only liq crosses bnd */
      if       (gasflowup==0&&gasflowdn==1) {
        /* unphysical case */
        if(jv*uval<0.0) {
          f = 0.0;
        } else {
          f = dSy(i,j,k) * jv * dt;
        }
        boil::aout<<"case1y "<<i<<" "<<jup<<" "<<k<<" | "<<jdn<<" "<<phiup<<" "<<stmp[i][jup][k]/dV(i,jup,k)<<" "<<f/dV(i,jup,k)<<" | "<<(stmp[i][jdn][k]+fabs(f))/dV(i,jdn,k)<<" "<<(stmp[i][jup][k]-fabs(f))/dV(i,jup,k)<<boil::endl;
      /* case 2: both downwind and upwind follow the flow */
      } else if(gasflowup==1||gasflowdn==1) { 
        /* calculate gliq and ggas: CFL of liquid and gas upwind */
        real gliq = uval*dt/dyup;
        real ggas;
        if(gasflowup==1) {
          ggas = gasvelup*dt/dyup;
        } else {
          ggas = gasveldn*dt/dyup;
        }
        /* dimensionless vol. flow = jv*dt*dS/dV*volfactor */
        real gj = jv*dt/dyup/volfactorup; 
#if 1
        f = dV(i,jup,k)*volfactorup
            *calc_diabatic_flux(gj,gliq,ggas,phiup,
                                ny[i][jup][k],nx[i][jup][k],nz[i][jup][k]);
        boil::aout<<"case2y "<<i<<" "<<jup<<" "<<k<<" | "<<jdn<<" "<<phiup<<" "<<stmp[i][jup][k]/dV(i,jup,k)<<" "<<f/dV(i,jup,k)<<" | "<<(stmp[i][jdn][k]+fabs(f))/dV(i,jdn,k)<<" "<<(stmp[i][jup][k]-fabs(f))/dV(i,jup,k)<<boil::endl;
#endif

      /* case 3: both downwind and upwind oppose the flow */ 
      } else if(gasflowdn==0||gasflowup==0) {
        real dydn = phi.dyc(jdn);

        /* degenerate case */
        if(dydn==0.0) {
          real jgas;
          if(gasflowdn==0) {
            jgas = (1.0-phi[i][jdn][k])*gasveldn;
          } else {
            jgas = (1.0-phi[i][jdn][k])*gasvelup;
          }
          jv -= jgas;

          /* unphysical case */
          if(jv*uval<0.0) {
            f = 0.0;
          } else {
            f = dSy(i,j,k)*jv*dt;
          } 
          boil::aout<<"case3y "<<i<<" "<<jup<<" "<<k<<" | "<<jdn<<" "<<phiup<<" "<<stmp[i][jup][k]/dV(i,jup,k)<<" "<<f/dV(i,jup,k)<<" | "<<(stmp[i][jdn][k]+fabs(f))/dV(i,jdn,k)<<" "<<(stmp[i][jup][k]-fabs(f))/dV(i,jup,k)<<boil::endl;
        } else {
          /* calculate gliq and ggas: CFL of liquid upw and gas dnw */
          real gliq = uval*dt/dyup;
          real ggas;
          if(gasflowdn==0) {
            ggas = gasveldn*dt/dydn;
          } else {
            ggas = gasvelup*dt/dydn;
          }
          real gj = jv*dt/dyup/volfactorup;

          /* calculate ratio of dy, including the effect of vol expansion */
          real dyrat = dydn/dyup*volfactordn/volfactorup;
#if 1
          f = volfactorup*dV(i,jup,k)
              *calc_diabatic_flux(gj,gliq,ggas,dyrat,phiup,phidn,
                                  ny[i][jup][k],nx[i][jup][k],nz[i][jup][k],
                                  ny[i][jdn][k],nx[i][jdn][k],nz[i][jdn][k]);
          boil::aout<<"case3y "<<i<<" "<<jup<<" "<<k<<" | "<<jdn<<" "<<phiup<<" "<<stmp[i][jup][k]/dV(i,jup,k)<<" "<<f/dV(i,jup,k)<<" | "<<(stmp[i][jdn][k]+fabs(f))/dV(i,jdn,k)<<" "<<(stmp[i][jup][k]-fabs(f))/dV(i,jup,k)<<boil::endl;
#endif
        }
      /* case 4: downwind opposes flow, upwind follows: all three cross bnd */
      } else {
        /* underdevelopment */
        boil::oout<<"VOF_advance_y: underdevelopment! Exiting."<<boil::endl;
        exit(0);
#if 0

#endif
      } /* cases */
    } /* adiabatic-diabatic */

#if 0
    /* limit flux */
    real sgnf = (f>0.0)-(f<0.0);
    f = sgnf*std::min(fabs(stmp[i][jup][k]),fabs(f));
#endif

#if 0
    /* update stmp */
    stmp[i][j-1][k] = stmp[i][j-1][k] - f;
    stmp[i][j  ][k] = stmp[i][j  ][k] + f;
#else
    fluxmax[m][i][j][k] = f;
    sosflux[m][i][j][k] = f/3.0;
#endif

#if 0
    if(i==100&&j==3&&k==100) {
      std::cout<<"advance_y:"<<f<<" "<<j<<" "<<phi[i][j][k]<<"\n";
    }
#endif

  }

  return;

}

