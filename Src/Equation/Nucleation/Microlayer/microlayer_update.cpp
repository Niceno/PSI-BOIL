#include "microlayer.h"

/******************************************************************************/
void Microlayer::update(real & smdot_micro,
                        real & smdot_pos_macro_overwrite,
                        real & smdot_neg_macro_overwrite) {
/***************************************************************************//**
*  \brief update microlayer thickness and heat source
*******************************************************************************/

#ifdef DEBUG
  std::cout<<"microlayer.update: "<<boil::cart.iam()<<"\n";
#endif

  const real dt = time->dt();

#ifndef USE_VOF
  /* reset heat sink */
  *tprs = 0.0;

  /* store area of vapor to dSprev, if not stored */
  if (!str_dSprev) {
    store_dSprev();
  }
#endif

  /* initialize sum variables */
  smdot_micro = 0.0;
  smdot_pos_macro_overwrite = 0.0;
  smdot_neg_macro_overwrite = 0.0;
  for(int i=0; i<=6; i++){
    area_sum[i]=0.0;
    area_l[i]=0.0;
    area_v[i]=0.0;
    hflux_micro[i]=0.0;
  }

  /*---------------------------------+
  |  calculate mdot on wall boundary |
  +---------------------------------*/
  for( int b=0; b<dmicro.bc().count(); b++ ) {
    if(dmicro.bc().type_decomp(b)) continue;
    if( dmicro.bc().type(b) == BndType::wall()) {
      int iof=0, jof=0, kof=0;
      Dir d = dmicro.bc().direction(b);
      int ndir = int(d);
      if(d != Dir::undefined()) {

        Comp mcomp;
        Sign sig;

        if (d == Dir::imin()) {
          mcomp = Comp::i();
          sig = Sign::neg();
          iof++;
        } else if (d == Dir::imax()) {
          mcomp = Comp::i();
          sig = Sign::pos();
          iof--;
        } else if (d == Dir::jmin()) {
          mcomp = Comp::j();
          sig = Sign::neg();
          jof++;
        } else if (d == Dir::jmax()) {
          mcomp = Comp::j();
          sig = Sign::pos();
          jof--;
        } else if (d == Dir::kmin()) {
          mcomp = Comp::k();
          sig = Sign::pos();
          kof++;
        } else if (d == Dir::kmax()) {
          mcomp = Comp::k();
          sig = Sign::pos();
          kof--;
        } else {
          continue;
        }

        for_vijk( dmicro.bc().at(b), i,j,k ){
          // skip solid cells
          if(cht->topo->domain()->ibody().off(i+iof,j+jof,k+kof)) continue;

          // only fluid cell comes here
          /*--------------------------------------------------+
          |  Note:  (i    , j    , k    ) is in solid domain  |
          |         (i+iof, j+jof, k+kof) is in fluid domain  |
          +--------------------------------------------------*/
          int ii=i+iof;
          int jj=j+jof;
          int kk=k+kof;
          real vol = dmicro.dV(ii,jj,kk);

          /* area */
          real area = std::abs(iof)*dmicro.dSx(sig,ii,jj,kk)
                    + std::abs(jof)*dmicro.dSy(sig,ii,jj,kk)
                    + std::abs(kof)*dmicro.dSz(sig,ii,jj,kk);
          area_sum[ndir] += area;

          real area_vap = area_vapor(sig,mcomp,ii,jj,kk);

          area_l[ndir] += area - area_vap;
          area_v[ndir] += area_vap;

#ifdef USE_VOF
          if(in_vapor(ii,jj,kk)) {
#else
          if(area_vap>0.) {
#endif
 
            /* micro layer model */
            if(!boil::realistic(dmicro[ii][jj][kk])) {
              dmicro[ii][jj][kk] = d0(ii,jj,kk);
            }

            /* micro layer does not exist = dry patch */
            if(dmicro[ii][jj][kk] <= dmicro_min*(1.+boil::pico)) {
	    } else {

              /* micro layer exists */
	      real dmicro_new;
              dmicro_new = dmicro[ii][jj][kk]
                         - dt/rhol/latent
                         * ((cht->tmp())[i][j][k] - cht->tifmodel.Tint(i,j,k)) 
                         / (dmicro[ii][jj][kk]/lambdal + hresis);
	      dmicro_new = std::max(dmicro_new, dmicro_min);

              /* micro-layer thickness doesn't increase */
              dmicro_new = std::min(dmicro_new,  dmicro[ii][jj][kk]); 

              real Mdot = -rhol * (dmicro_new - dmicro[ii][jj][kk])
                        / dt * area / vol;
              (*mdot)[ii][jj][kk] = Mdot;
	      dmicro[ii][jj][kk]=dmicro_new;

              /* heat flux */
              hflux_micro[ndir] += Mdot * vol/area * latent;

            }

          } /* in vapor */
        } /* all cells */
      } /* dir not undefined */
    } /* is wall */
  } /* all bcs */

  /*-------------------------------------+
  |  calculate mdot on immersed boundary |
  +-------------------------------------*/
  for(int cc=0; cc<cht->topo->domain()->ibody().nccells(); cc++){
    int i,j,k;
    cht->topo->domain()->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=cht->topo->domain()->ibody().nwx(i,j,k);
    real uy=cht->topo->domain()->ibody().nwy(i,j,k);
    real uz=cht->topo->domain()->ibody().nwz(i,j,k);
    Dir d = Dir::undefined();
    Sign sig;
    Comp mcomp;
    if (fabs(uz)>0.707) {
      mcomp = Comp::k();
      d = Dir::kmin();
      if(uz>0.707){
        d = Dir::kmax();
        sig = Sign::pos();
      }
    } else if (fabs(ux)>0.707) {
      mcomp = Comp::i();
      d = Dir::imin();
      if(ux>0.707){
        d = Dir::imax();
        sig = Sign::pos();
      }
    } else if (fabs(uy)>0.707) {
      mcomp = Comp::j();
      d = Dir::jmin();
      if(uy>0.707){
        d = Dir::jmax();
        sig = Sign::pos();
      }
    } else {
      boil::oout<<"microlayer_update: Underdevelopment!!!\n";
      exit(0);
    }

    int iof=0, jof=0, kof=0;
    int iofv=0, jofv=0, kofv=0;
    int ndir = 6;  // ibody()
    if(d == Dir::imin()) { iof--; }; 
    if(d == Dir::imax()) { iof++; iofv++; };
    if(d == Dir::jmin()) { jof--; }; 
    if(d == Dir::jmax()) { jof++; jofv++; };
    if(d == Dir::kmin()) { kof--; }; 
    if(d == Dir::kmax()) { kof++; kofv++; };

    /*--------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in fluid domain  |
    |         (i+iof, j+jof, k+kof) is in solid domain  |
    +--------------------------------------------------*/
    real vol = dmicro.dV(i,j,k);

    /* area */
    real area = std::abs(iof)*dmicro.dSx(sig,i,j,k)
              + std::abs(jof)*dmicro.dSy(sig,i,j,k)
              + std::abs(kof)*dmicro.dSz(sig,i,j,k);
    area_sum[ndir] += area;

    real area_vap = area_vapor(sig,mcomp,i,j,k);
    area_l[ndir] += area - area_vap;
    area_v[ndir] += area_vap;

    int ii=i+iof;
    int jj=j+jof;
    int kk=k+kof;

#ifdef USE_VOF
    if(!in_vapor(i,j,k)) {
#else
    if(area_vap == 0.0) {
#endif
    } else {
#ifndef USE_VOF
      /* here, the area of the microlayer is not considered!!! 
       * this would require correction */
      exit(0);
      //if ( approx (area_vap, area, area*boil::micro)) 
#endif

      if ( dmicro[i][j][k] <= dmicro_min*(1.+boil::pico) // depleted
        || !boil::realistic(dmicro[i][j][k])) {          // full vapor
      } else {

        /* micro layer exists */
        real dmicro_new;
        real dt = time->dt();

        real qtmp = ((cht->node_tmp_flu())[mcomp][i+iofv][j+jofv][k+kofv]
                    - cht->tifmodel.Tint(i,j,k))
                  / ( dmicro[i][j][k]/lambdal + hresis ); 

        dmicro_new = dmicro[i][j][k] - dt / rhol * qtmp / latent;
        dmicro_new = std::max(dmicro_new, dmicro_min);

        /* micro-layer thickness doesn't increase */
        dmicro_new = std::min(dmicro_new, dmicro[i][j][k]); 

        if ((*mdot)[i][j][k]>=0.0) {
          smdot_pos_macro_overwrite += (*mdot)[i][j][k]*vol;
        } else {
          smdot_neg_macro_overwrite += (*mdot)[i][j][k]*vol;
        }

        /* overwrite mdot */
#ifdef USE_VOF
        real mdot_micro = - rhol * (dmicro_new - dmicro[i][j][k])
                          / dt * area / vol;
        (*mdot)[i][j][k] = mdot_micro;
#else
        real mdot_micro = - rhol * (dmicro_new - dmicro[i][j][k])
                          / dt * area_vap / vol;
        (*mdot)[i][j][k] += mdot_micro;
#endif
        smdot_micro += mdot_micro*vol;

#ifndef USE_VOF
        /* removed due to update-at-walls */
        /* enthalpy clean-up: sink due to microlayer */
        (*tprs)[i+iof][j+jof][k+kof] += -vol*mdot_micro*latent;
        if(in_vapor(i,j,k)) {
          /* additional effect due to heat-up in vapour */
          (*tprs)[i][j][k] -= cpv *((cht->tmp())[i][j][k]-cht->tifmodel.Tint(i,j,k))
                              //*(1.0/rhov-1.0/rhol)*mdot_micro*vol;
                              *1.0/rhov*mdot_micro*vol;
        }
#endif

        /* microlayer thickness changes */
        dmicro[i][j][k]=dmicro_new;

        /* heat flux */
        /* it doesn't take into account max(dmicro_new, dmicro_min) */
        hflux_micro[ndir] += qtmp*area_vap;
        /* it does take into account max(dmicro_new, dmicro_min) */
        //hflux_micro[ndir] += (*mdot)[i][j][k] * vol * latent;

      } /* not zero microlayer thickness */
    } /* microlayer possible */

  } /* ibody cells */

  mdot->exchange_all();

  boil::cart.sum_real_n(area_sum,7);
  boil::cart.sum_real_n(area_l,7);
  boil::cart.sum_real_n(area_v,7);
  boil::cart.sum_real_n(hflux_micro,7);

#ifdef DEBUG
  boil::plot->plot(dmicro, (cht->tmp()), *mdot, "dmicro-tpr-mdot",  time->current_step());
#endif

  boil::cart.sum_real(&smdot_micro);
  boil::cart.sum_real(&smdot_pos_macro_overwrite);
  boil::cart.sum_real(&smdot_neg_macro_overwrite);

  boil::oout<<"micro_overwrite: time= "<<time->current_time()
            <<" smdot_micro= "<<smdot_micro
            <<" smdot_pos_overwrite= "<<smdot_pos_macro_overwrite
            <<" smdot_neg_overwrite= "<<smdot_neg_macro_overwrite
            <<"\n";


#ifdef DEBUG
  boil::plot->plot(*cht->topo->clr, dmicro, *mdot, 
                   "clr-dmicro-mdot",  time->current_step());
  exit(0);
#endif

#ifndef USE_VOF
  /* store area of vapor to dSprev */
  store_dSprev();
#endif

  return;
}
