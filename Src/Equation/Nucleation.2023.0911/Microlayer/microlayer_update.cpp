#include "microlayer.h"
#include "../header.h"

/******************************************************************************/
void Microlayer::update(const Scalar * diff_eddy) {
  real r1,r2,r3;
  update(r1,r2,r3,diff_eddy,NULL,NULL);
}
/******************************************************************************/
void Microlayer::update(Scalar * hflux, Scalar * warea){
  real r1,r2,r3;
  update(r1,r2,r3,NULL,hflux,warea);
}
/******************************************************************************/
void Microlayer::update(const Scalar * diff_eddy,
                        Scalar * hflux, Scalar * warea){
  real r1,r2,r3;
  update(r1,r2,r3,diff_eddy,hflux,warea);
}
/******************************************************************************/
void Microlayer::update(real & smdot_micro, const Scalar * diff_eddy) {
  real r2,r3;
  update(smdot_micro,r2,r3,diff_eddy,NULL,NULL);
}
/******************************************************************************/
void Microlayer::update(real & smdot_micro,
                        Scalar * hflux, Scalar * warea){
  real r2,r3;
  update(smdot_micro,r2,r3,NULL,hflux,warea);
}
/******************************************************************************/
void Microlayer::update(real & smdot_micro, const Scalar * diff_eddy,
                        Scalar * hflux, Scalar * warea){
  real r2,r3;
  update(smdot_micro,r2,r3,diff_eddy,hflux,warea);
}
/******************************************************************************/
void Microlayer::update(real & smdot_micro,
                        real & smdot_pos_macro_overwrite,
                        real & smdot_neg_macro_overwrite,
                        const Scalar * diff_eddy){
  real r1,r2,r3;
  update(smdot_micro, smdot_pos_macro_overwrite,smdot_neg_macro_overwrite,
         diff_eddy,NULL,NULL);
}
/******************************************************************************/
void Microlayer::update(real & smdot_micro,
                        real & smdot_pos_macro_overwrite,
                        real & smdot_neg_macro_overwrite,
                        Scalar * hflux,
                        Scalar * warea){
  real r1,r2,r3;
  update(smdot_micro, smdot_pos_macro_overwrite,smdot_neg_macro_overwrite,
         NULL,hflux,warea);
}

/******************************************************************************/
void Microlayer::update(real & smdot_micro,
                        real & smdot_pos_macro_overwrite,
                        real & smdot_neg_macro_overwrite,
                        const Scalar * diff_eddy,
                        Scalar * hflux,
                        Scalar * warea) {
/***************************************************************************//**
*  \brief update microlayer thickness and heat source
*
*  smdot_micro: sum of mdot in micro-layer [kg/s]
*  smdot_pos_macro_overwrite: mdot overwirtten in this function (positive part)
*  smdot_neg_macro_overwrite: mdot overwritten in this function (negative part)
*  diff_eddy: eddy viscosity
*  hflux: store heat flux at wall
*  warea: store area at wall  0~1: ratio of liquid area
*                            -1~0: -(ratio of micro-layer)
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"microlayer_update: "<<boil::cart.iam()<<"\n";
  if (diff_eddy) boil::oout<<"microlayer_update:diff_eddy\n";
  if (hflux) boil::oout<<"microlayer_update:hflux\n";
  if (warea) boil::oout<<"microlayer_update:warea\n";
#endif

  const real dt = time->dt();

#ifndef USE_VOF_NUCL
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
    area_micro[i]=0.0;
    hflux_sum[i]=0.0;
    hflux_l[i]=0.0;
    hflux_v[i]=0.0;
    hflux_micro[i]=0.0; 
  }

  // In the original version (Github), dmicro_intermediate() is called here.
  // But it is neglected here, because it must be trivial.

  /*----------------+
  |  wall boundary  |
  +----------------*/
  for( int b=0; b<dmicro.bc().count(); b++ ) {
    if(dmicro.bc().type_decomp(b)) continue;
    if( dmicro.bc().type(b) == BndType::wall()) {
      boil::oout<<"microlayer_update: underconstruction wall boundary\n";
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
          real area = std::abs(iof)*dmicro.dSx(sig,ii,jj,kk)
                    + std::abs(jof)*dmicro.dSy(sig,ii,jj,kk)
                    + std::abs(kof)*dmicro.dSz(sig,ii,jj,kk);
          real area_vap = area_vapor(sig,mcomp,ii,jj,kk);

          /*---------------------+
          |  statistics of area  |
	  +---------------------*/
          // total
          area_sum[ndir] += area;

          // liquid
          area_l[ndir]   += area - area_vap;
          if (warea) (*warea)[ii][jj][kk] = (area-area_vap)/area;

          // vapor
          if (dmicro[ii][jj][kk] >boil::mega) {
            // liquid + vapor (no micro-layer in vapor) -> no micro-layer
            area_v[ndir] += area_vap;
          } else if (dmicro[ii][jj][kk] <= dmicro_min*(1.+boil::pico)) {
            // liquid + vapor (dry spot) -> no micro-layer
            area_v[ndir] += area_vap;
          } else {
            // liquid + vapor (with microlayer)
            area_micro[ndir] += area_vap;
            if (warea) (*warea)[ii][jj][kk] = -area_vap/area; // -(ratio micro)
          }

          /*-----------------------------------+
          |  calculate dmicro, mdot and hflux  |
	  +-----------------------------------*/
#ifdef USE_VOF_NUCL
          if(below_threshold(ii,jj,kk)) {
#else
          if(area_vap == 0.0) {
            /* 100% liquid cell */
            // dmicro
            dmicro[ii][jj][kk]=boil::unreal;
          } else if (area_vap>0.0) {
#endif
            /* vapor exists in the cell */

            /* create new micro layer model */
            if(!boil::realistic(dmicro[ii][jj][kk])) { //(=dmicro > mega)
              /* If micro-layer didn't exist, create new micro-layer */
              dmicro[ii][jj][kk] = d0(ii,jj,kk);
            } else if(dmicro[ii][jj][kk] <= dmicro_min*(1.+boil::pico)) {
              /* dry-patch */
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
#ifndef USE_VOF_NUCL
              real mdot_micro = -rhol * (dmicro_new - dmicro[ii][jj][kk])
                                / dt * area / vol;
#else
              real mdot_micro = -rhol * (dmicro_new - dmicro[ii][jj][kk])
                                / dt * area_vap / vol;
#endif
              (*mdot)[ii][jj][kk] = mdot_micro;
              smdot_micro += mdot_micro*vol;
	      dmicro[ii][jj][kk]=dmicro_new;

            }
          }
        } /* all cells */
      } /* dir not undefined */
    } /* is wall */
  } /* all bcs */

  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
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
    Sign sig = Sign::neg();
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
    int ii=i+iof;
    int jj=j+jof;
    int kk=k+kof;
    real alen = fabs(dmicro.xc(i)-dmicro.xc(ii))
              + fabs(dmicro.yc(j)-dmicro.yc(jj))
              + fabs(dmicro.zc(k)-dmicro.zc(kk));
    // clude code
    real dw = 0.5 * dmicro.dzc(kk);  // half cell size in wall
    real df = 0.5 * dmicro.dzc(k);   // half cell size in fluid
    real lambdas = cht->solid()->lambda(ii,jj,kk);

    real vol = dmicro.dV(i,j,k);
    real area = std::abs(iof)*dmicro.dSx(sig,i,j,k)
              + std::abs(jof)*dmicro.dSy(sig,i,j,k)
              + std::abs(kof)*dmicro.dSz(sig,i,j,k);
    real area_vap = area_vapor(sig,mcomp,i,j,k);
    real area_liq = area - area_vap;

    /*---------------------+
    |  statistics of area  |
    +---------------------*/
    // total
    area_sum[ndir] += area;

    // liquid
    area_l[ndir]   += area_liq;
    if (warea) (*warea)[ii][jj][kk] = area_liq/area;

    // vapor
    if (dmicro[i][j][k] >boil::mega) {
      // liquid + vapor (no micro-layer in vapor) -> no micro-layer
      if(area_vap>0)
      //std::cout<<"micro_update:mega= "<<dmicro[i][j][k]<<" "<<area_vap<<" "
      //   <<area<<" "<<ii<<" "<<jj<<" "<<kk<<" "<<dmicro.zc(kk)<< "\n";
      area_v[ndir] += area_vap;
    } else if (dmicro[i][j][k] <= dmicro_min*(1.+boil::pico)) {
      //std::cout<<"micro_update:dry= "<<dmicro[i][j][k]<<" "<<dmicro_min<<"\n";
      // liquid + vapor (dry spot) -> no micro-layer
      area_v[ndir] += area_vap;
    } else {
      // liquid + vapor (with microlayer)
      //std::cout<<"micro_update:micro= "<<dmicro[i][j][k]<<"\n";
      area_micro[ndir] += area_vap;
      if (warea) (*warea)[ii][jj][kk] = -area_vap/area; // -(ratio micro)
    }


    /*-----------------------------------+
    |  calculate dmicro, mdot and hflux  |
    +-----------------------------------*/
#ifdef USE_VOF_NUCL
    //if(!in_vapor(i,j,k)) {
    if(!below_threshold(i,j,k)) {
    } else {
#else
    if(area_vap == 0.0) {
      /*-------------------+
      |  100% liquid cell  |
      +-------------------*/

      /* dmicro calculation */
      dmicro[i][j][k]=boil::unreal;

      /* hflux calculation */
      real laml = cht->lambdal(i,j,k,diff_eddy);
      real twl = ( df * lambdas * cht->tmp()[ii][jj][kk]
                 + dw * laml    * cht->tmp()[i ][j ][k ] )
                 / ( df * lambdas + dw * laml);
      real ql = lambdas*(cht->tmp()[ii][jj][kk]-twl)/dw; //[W/m2]
      hflux_sum[ndir] += ql*area; //[W]
      hflux_l[ndir] += ql*area; //[W]
      if(hflux) (*hflux)[ii][jj][kk]=ql; //[W/m2]

    } else if (area_vap>0.0) {
#endif
      /*---------------------------------------+
      |  100% vapor cell or liquid-vapor cell  |
      +---------------------------------------*/
      if ( dmicro[i][j][k] <= dmicro_min*(1.+boil::pico) // depleted
        || !boil::realistic(dmicro[i][j][k])) {          // no micro-layer
        /*-----------------+
        |  no micro-layer  |
        +-----------------*/
        /* hflux calculation */
	// through vapor phase
        real lamv = cht->lambdav(i,j,k,diff_eddy);
        real twv = ( df * lambdas * cht->tmp()[ii][jj][kk]
                   + dw * lamv    * cht->tmp()[i ][j ][k ] )
                  / ( df * lambdas + dw * lamv);
        real qv = lambdas*(cht->tmp()[ii][jj][kk]-twv)/dw; //[W/m2]
	// through liquid phase (if it exists)
        real laml = cht->lambdal(i,j,k,diff_eddy);
        real twl = ( df * lambdas * cht->tmp()[ii][jj][kk]
                   + dw * laml    * cht->tmp()[i ][j ][k ] )
                 / ( df * lambdas + dw * laml);
        real ql = lambdas*(cht->tmp()[ii][jj][kk]-twl)/dw; //[W/m2]

        hflux_sum[ndir] += qv*area_vap + ql*area_liq; // [W]
        hflux_l[ndir] += ql*area_liq; // [W]
        hflux_v[ndir] += qv*area_vap; // [W]
        if (hflux) 
          (*hflux)[ii][jj][kk] = (qv*area_vap+ql*area_liq)/area;//[W/m2]

      } else {
        /*---------------------+
        |  micro layer exists  |
        +---------------------*/
        /* dmicro calculation */
#ifdef USE_VOF_NUCL
        real qmic = ((cht->node_tmp_flu())[mcomp][i+iofv][j+jofv][k+kofv]
                    - cht->tifmodel.Tint(i,j,k))
                  / ( dmicro[i][j][k]/lambdal + hresis ); // hresis: resistance
#else
        real qmic = (cht->tmp()[ii][jj][kk] - cht->tifmodel.Tint(i,j,k))
                  / ( dw/lambdas + dmicro[i][j][k]/lambdal + hresis);
#endif
	if (qmic<0) {
          boil::aout<<"### WARNING:microlayer_update:negative_qmic= "<<qmic<<" proc= "
          <<boil::cart.iam()<<" i " <<i<<" j "<<j<<" k "<<k<<" kk "<<kk<<" T_interface "
          <<cht->tmp()[ii][jj][kk]<<" Tw "<<cht->tifmodel.Tint(i,j,k)<<" dmicro "
          <<dmicro[i][j][k]<<"\n";
          //exit(0);
	}

        real dmicro_new;
        dmicro_new = dmicro[i][j][k] - dt / rhol * qmic / latent;
        dmicro_new = std::max(dmicro_new, dmicro_min);
        // micro-layer thickness doesn't increase
        dmicro_new = std::min(dmicro_new, dmicro[i][j][k]);

        // change of micro-layer thickness in this time step
	real delta_micro = dmicro[i][j][k] - dmicro_new;

        // statistics before overwritten
        if ((*mdot)[i][j][k]>=0.0) {
          smdot_pos_macro_overwrite += (*mdot)[i][j][k]*vol;
        } else {
          smdot_neg_macro_overwrite += (*mdot)[i][j][k]*vol;
        }

        /* mdot calculation */
#ifdef USE_VOF_NUCL
        real mdot_micro = - rhol * (dmicro_new - dmicro[i][j][k])
                          / dt * area / vol;
#else
        real mdot_micro = - rhol * (dmicro_new - dmicro[i][j][k])
                          / dt * area_vap / vol;
        //std::cout<<"microlayer_update:mdot= "<<i<<" "<<j<<" "<<k<<" "<<dmicro[i][j][k]
	//	<<" "<<dmicro_new<<" "<<qmic<<" "<<latent<<"\n";
#endif
        (*mdot)[i][j][k] = mdot_micro;
        smdot_micro += mdot_micro*vol;

        // microlayer thickness changes
        dmicro[i][j][k]=dmicro_new;

#ifndef USE_VOF_NUCL
        /* tprs calculation */
	// NEED TO CONSIDER 2023.01.18
        /* removed due to update-at-walls */
        /* enthalpy clean-up: sink due to microlayer */
        (*tprs)[i+iof][j+jof][k+kof] += -vol*mdot_micro*latent; // [W]
#if 0
       if(boil::cart.iam()==8&&i==21&&j==77){
          std::cout<<"micro_update:area_vap= "<<area_vap<<" "<<dmicro[i][j][k]<<" "
		  <<!boil::realistic(dmicro[i][j][k])<<"\n";
          std::cout<<"micro_update:dmicro= "<<dmicro[21][77][26]<<" "<<dmicro[21][77][27]<<" "
		   <<(*tprs)[i+iof][j+jof][k+kof]<<" "<<k<<"\n";
       }
#endif
        if(in_vapor(i,j,k)) {
          /* additional effect due to heat-up in vapour */
          (*tprs)[i][j][k] -= cpv * 
                              ((cht->tmp())[i][j][k]-cht->tifmodel.Tint(i,j,k))
                              *1.0/rhov*mdot_micro*vol;
                              //*(1.0/rhov-1.0/rhol)*mdot_micro*vol;
        }
#endif
        /* heat flux calculation */
        // through liquid phase (if it exists)
        real laml = cht->lambdal(i,j,k,diff_eddy);
        real twl = ( df * lambdas * cht->tmp()[ii][jj][kk]
                   + dw * laml    * cht->tmp()[i ][j ][k ] )
                 / ( df * lambdas + dw * laml);
        real ql = lambdas*(cht->tmp()[ii][jj][kk]-twl)/dw; //[W/m2]

        /* it doesn't take into account max(dmicro_new, dmicro_min) */
        hflux_micro[ndir] += qmic*area_vap;  // NEED TO CONSIDER 2023.01.18
        hflux_l[ndir] += ql*area_liq;
        hflux_sum[ndir] += qmic*area_vap + ql*area_liq;
        if (hflux)
          (*hflux)[ii][jj][kk] = (qmic*area_vap+ql*area_liq)/area;//[W/m2]

      }
    }
  } /* ibody cells */

  mdot->exchange_all();

  boil::cart.sum_real_n(area_sum,7);
  boil::cart.sum_real_n(area_l,7);
  boil::cart.sum_real_n(area_v,7);
  boil::cart.sum_real_n(area_micro,7);
  boil::cart.sum_real_n(hflux_sum,7);
  boil::cart.sum_real_n(hflux_l,7);
  boil::cart.sum_real_n(hflux_v,7);
  boil::cart.sum_real_n(hflux_micro,7);

  /* console-output */
  boil::oout<<"micro_update:area:[m2] "<<time->current_time()
            <<" total "<<area_sum[6]
            <<" lquid "<<area_l[6]
            <<" vapor "<<area_v[6]
            <<" micro "<<area_micro[6]<<"\n";
  boil::oout<<"micro_update:hflux:[W] "<<time->current_time()
            <<" total "<<hflux_sum[6]
            <<" lquid "<<hflux_l[6]
            <<" vapor "<<hflux_v[6]
            <<" micro "<<hflux_micro[6]<<"\n";

#ifdef DEBUG
  boil::plot->plot(dmicro, (cht->tmp()), *mdot, "dmicro-tpr-mdot",  time->current_step());
#endif

  boil::cart.sum_real(&smdot_micro);
  boil::cart.sum_real(&smdot_pos_macro_overwrite);
  boil::cart.sum_real(&smdot_neg_macro_overwrite);

#if 0
  boil::oout<<"micro_update:smdot: "<<time->current_time()
            <<" smdot_micro "<<smdot_micro
            <<" smdot_pos_overwrite "<<smdot_pos_macro_overwrite
            <<" smdot_neg_overwrite "<<smdot_neg_macro_overwrite<<"\n";
#endif

#ifdef DEBUG
  boil::plot->plot(*cht->topo->clr, dmicro, *mdot, 
                   "clr-dmicro-mdot",  time->current_step());
  exit(0);
#endif

#ifndef USE_VOF_NUCL
  /* store area of vapor to dSprev */
  store_dSprev();
#endif

  return;
}
