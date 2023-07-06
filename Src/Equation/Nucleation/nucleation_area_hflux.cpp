#include "nucleation.h"
#include "header.h"

/******************************************************************************/
void Nucleation::area_hflux(const Scalar * diff_eddy) {
  area_hflux(diff_eddy, NULL, NULL);
}
/******************************************************************************/
void Nucleation::area_hflux(Scalar * hflux, Scalar * warea){
  area_hflux(NULL, hflux, warea);
}

/******************************************************************************/
void Nucleation::area_hflux(const Scalar * diff_eddy,
                    Scalar * hflux, Scalar * warea) {
/***************************************************************************//**
*  \brief calculate area and hflux at boundaries
*  This kind of calculations should be done in CommonHeatTransfer, but it is 
*  done here, because of area_vapor.
*  area_vapor is the area of vapor computed with 2D marching cube.
*  CommonHeatTransfer doesn't have such a function.
*
*  If microlayer is defined, it is computed in micro.update();
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"nucleation_area_hflux: "<<boil::cart.iam()<<"\n";
  if (diff_eddy) boil::oout<<"microlayer_update:diff_eddy\n";
  if (hflux) boil::oout<<"microlayer_update:hflux\n";
  if (warea) boil::oout<<"microlayer_update:warea\n";
#endif

  const real dt = time->dt();

  /* initialize sum variables */
  real smdot_micro = 0.0;
  real smdot_pos_macro_overwrite = 0.0;
  real smdot_neg_macro_overwrite = 0.0;

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

  /*----------------+
  |  wall boundary  |
  +----------------*/
  for( int b=0; b<cht->tpr.bc().count(); b++ ) {
    if( cht->tpr.bc().type_decomp(b) ) continue;
    if( cht->tpr.bc().type(b) == BndType::wall()) {
      boil::oout<<"microlayer_update: underconstruction wall boundary\n";
      int iof=0, jof=0, kof=0;
      Dir d = cht->tpr.bc().direction(b);
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

        for_vijk( cht->tpr.bc().at(b), i,j,k ){
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
          real vol = cht->tpr.dV(ii,jj,kk);
          real area = std::abs(iof)*cht->tpr.dSx(sig,ii,jj,kk)
                    + std::abs(jof)*cht->tpr.dSy(sig,ii,jj,kk)
                    + std::abs(kof)*cht->tpr.dSz(sig,ii,jj,kk);
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
          area_v[ndir] += area_vap;

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
    real alen = fabs(cht->tpr.xc(i)-cht->tpr.xc(ii))
              + fabs(cht->tpr.yc(j)-cht->tpr.yc(jj))
              + fabs(cht->tpr.zc(k)-cht->tpr.zc(kk));
    // clude code
    real dw = 0.5 * cht->tpr.dzc(kk);  // half cell size in wall
    real df = 0.5 * cht->tpr.dzc(k);   // half cell size in fluid
    real lambdas = cht->solid()->lambda(ii,jj,kk);

    real vol = cht->tpr.dV(i,j,k);
    real area = std::abs(iof)*cht->tpr.dSx(sig,i,j,k)
              + std::abs(jof)*cht->tpr.dSy(sig,i,j,k)
              + std::abs(kof)*cht->tpr.dSz(sig,i,j,k);
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
    area_v[ndir] += area_vap;

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
    }
  } /* ibody cells */

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

  return;
}
