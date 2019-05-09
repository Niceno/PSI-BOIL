#include "vof.h"

/******************************************************************************/
real VOF::vel_value(const Comp m, const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief Calculate liquid velocity at given position.
    fext is related to mflx as follows:
    fext = -m'''/rhol = -m''*adens/rhol 
*******************************************************************************/
#if 1
   real uval = (*u)[m][i][j][k];
#else
   real uval = uliq[m][i][j][k];
#endif

#if 0
   if(mixture) {
     /* simplified approach */
     int ii(i), jj(j), kk(k);
     bool dirx(false), diry(false), dirz(false);
     if       (m == Comp::i()) {
       ii--;
       dirx = true;
     } else if(m == Comp::j()) {
       jj--;
       diry = true;
     } else {
       kk--;
       dirz = true;
     }
  
     real mflxm = fext[ii][jj][kk];
     real mflxp = fext[i ][j ][k ];

     real bdphi;
     if(bndclr) {
       bdphi = (*bndclr)[m][i][j][k];
     } else {
       bdphi = 0.5*(phi[ii][jj][kk]+phi[i][j][k]);
     }

     bool pcmin = fabs(mflxm)>boil::pico;
     bool pcplu = fabs(mflxp)>boil::pico;
 
  #if 1 /* for 1D receding film, turn this off */
     /* phase change occurs in both cells: uval assumed to be the vol. avg */
     if(pcmin&pcplu) {
     #if 1
       real coef = (1.0-bdphi)*(rhol/rhov-1.0); 
       uval += vel_correct(ii,jj,kk,dirx,diry,dirz,coef,mflxm);
       uval += vel_correct(i ,j ,k ,dirx,diry,dirz,coef,mflxp);

       //boil::aout<<"M+P: "<<m<<" "<<ii<<" "<<jj<<" "<<kk<<" "<<mflxm<<" "<<bdphi<<" "<<vel_correct(ii,jj,kk,dirx,diry,dirz,coef,mflxm)<<boil::endl;
       //boil::aout<<"M+P: "<<m<<" "<<i<<" "<<j<<" "<<k<<" "<<mflxp<<" "<<bdphi<<" "<<vel_correct(i,j,k,dirx,diry,dirz,coef,mflxp)<<boil::endl;
     #endif
     /* phase change occurs in the minus cell: uval assumed to be single phase */ 
     } else if(pcmin) {
       /* staggered cell center does not correspond to normal cell boundary */
       //if(bdphi<phisurf) { 
       if(phi[i ][j ][k ]<phisurf) {
         real coef = (rhol/rhov-1.0);
         uval += vel_correct(ii,jj,kk,dirx,diry,dirz,coef,mflxm);
         //boil::aout<<"M: "<<m<<" "<<ii<<" "<<jj<<" "<<kk<<" "<<mflxm<<" "<<bdphi<<" "<<phi[i][j][k]<<" "<<vel_correct(ii,jj,kk,dirx,diry,dirz,coef,mflxm)<<boil::endl;
       }
     /* phase change occurs in the plus cell: uval assumed to be single phase */ 
     } else if(pcplu) {
       //if(bdphi<phisurf) {
       if(phi[ii][jj][kk]<phisurf) {
         real coef = (rhol/rhov-1.0);
         uval += vel_correct(i ,j ,k ,dirx,diry,dirz,coef,mflxp);
         //boil::aout<<"P: "<<m<<" "<<i<<" "<<j<<" "<<k<<" "<<mflxp<<" "<<bdphi<<" "<<phi[ii][jj][kk]<<" "<<vel_correct(i,j,k,dirx,diry,dirz,coef,mflxp)<<boil::endl;
       }
     }
  #else
     if(pcmin) {
       /* staggered cell center does not correspond to normal cell boundary */
       //if(bdphi<phisurf) { 
       if(phi[i ][j ][k ]<phisurf) {
         real coef = (rhol/rhov-1.0);
         uval += vel_correct(ii,jj,kk,dirx,diry,dirz,coef,mflxm);
         //boil::aout<<"M: "<<m<<" "<<ii<<" "<<jj<<" "<<kk<<" "<<mflxm<<" "<<bdphi<<" "<<phi[i][j][k]<<" "<<vel_correct(ii,jj,kk,dirx,diry,dirz,coef,mflxm)<<boil::endl;
       }
     } else if(pcplu) {
       //if(bdphi<phisurf) {
       if(phi[ii][jj][kk]<phisurf) {
         real coef = (rhol/rhov-1.0);
         uval += vel_correct(i ,j ,k ,dirx,diry,dirz,coef,mflxp);
         //boil::aout<<"P: "<<m<<" "<<i<<" "<<j<<" "<<k<<" "<<mflxp<<" "<<bdphi<<" "<<phi[ii][jj][kk]<<" "<<vel_correct(i,j,k,dirx,diry,dirz,coef,mflxp)<<boil::endl;
       }
     }
  #endif
   } /* if mixture */
#endif

   return uval;
}

real VOF::vel_correct(const int i, const int j, const int k,
                      const bool dirx, const bool diry, const bool dirz,
                      const real coef, const real mflx) {
  real mmx = nx[i][j][k];
  real mmy = ny[i][j][k];
  real mmz = nz[i][j][k];

  real nnx = mmx/phi.dxc(i);
  real nny = mmy/phi.dyc(j);
  real nnz = mmz/phi.dzc(k);

  real nnorm = nnx*nnx+nny*nny+nnz*nnz;
  nnorm = sqrt(nnorm)+boil::pico;

  /* aligned direction */
  real ndir = dirx*nnx+diry*nny+dirz*nnz;
  ndir /= nnorm;
#if 0
  boil::oout<<i<<" "<<j<<" "<<k<<" "<<-coef*mflx*ndir/(adens[i][j][k]+boil::pico)<<boil::endl;
#endif
  /* - due to the direction of the normal vector */
  return -coef*mflx*ndir/(adens[i][j][k]+boil::pico);
}
