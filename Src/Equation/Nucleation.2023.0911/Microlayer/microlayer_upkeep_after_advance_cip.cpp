#include "microlayer.h"
#ifndef USE_VOF_NUCL

/******************************************************************************/
void Microlayer::upkeep_after_advance_cip() {
/***************************************************************************//**
*  \brief update microlayer thickness due to effect of bubble area change
*******************************************************************************/
  //boil::oout<<"upkeep_cip: begin\n";

  for(int cc=0; cc<cht->topo->domain()->ibody().nccells(); cc++){
    int i,j,k;
    cht->topo->domain()->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=cht->topo->domain()->ibody().nwx(i,j,k);
    real uy=cht->topo->domain()->ibody().nwy(i,j,k);
    real uz=cht->topo->domain()->ibody().nwz(i,j,k);
    Comp mcomp;
    Sign sig = Sign::neg();
    Dir d = Dir::undefined();
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
      boil::oout<<"microlayer_area_effect: Underdevelopment!!!\n";
      exit(0);
    }

    int iof=0, jof=0, kof=0;
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*--------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in fluid domain  |
    |         (i+iof, j+jof, k+kof) is in solid domain  |
    +--------------------------------------------------*/
    int ii=i+iof;
    int jj=j+jof;
    int kk=k+kof;

    /* include vapor */
    real a_vapor=area_vapor(sig,mcomp,i,j,k);
    //boil::oout<<"upkeep_cip: "<<i<<" "<<j<<" "<<k<<" "<<a_vapor<<"\n";
#if 0
    if (i==9&&j==10&&k==15) {
      std::cout<<"upkeep_cip:9,10,15 "<<a_vapor<<" "
               <<(*cht->topo->clr)[i][j][k]<<"\n";
    }
#endif
    //if (a_vapor>0.0) std::cout<<"upkeep_cip:a_vapor= "<<a_vapor
    //                          <<" "<<i<<" "<<j<<" "<<k<<"\n";

    if(a_vapor > 0.0) {

      /* initialize */
      if(!boil::realistic(dmicro[i][j][k])) {
        dmicro[i][j][k] = d0(i,j,k);
        // vf = dmicro *dSz / vol = dmicro / dzc 
#if 0
        if(i==23 && j==23 && k==15){
          std::cout<<"upkeep_cip:23,23,15 "<<dmicro[i][j][k]<<" "
		  <<cht->distance_face(sig,mcomp,i,j,k)<<"\n";
	}
#endif
#if 0
        (*vf)[i][j][k] = dmicro[i][j][k]/(2.*cht->distance_face(sig,mcomp,i,j,k));
        (*cht->topo->clr)[i][j][k] = dmicro[i][j][k]/(2.*cht->distance_face(sig,mcomp,i,j,k));
#endif

      }

      real area = dmicro.dSz(sig,i,j,k);

      /* finalize to prevent microlayer regeneration */
      if( (dmicro[i][j][k]<=dmicro_min*(1.+boil::pico)) &&
          (approx(a_vapor, area, area*boil::micro))
        ) {
        dmicro[i][j][k] = 0.0;
      }

      /* interface cell */
      if(dSprev[i][j][k] < area) {
        real s_prev = dSprev[i][j][k];
        real s_now  = a_vapor;
        if( (s_now > s_prev) &&
            //(dmicro[i][j][k] > dmicro_min*(1.+boil::pico)) &&
            (dmicro[i][j][k] > 0.0) &&
            boil::realistic(dmicro[i][j][k])
          ) {
          real d_new = ( dmicro[i][j][k] * s_prev
                       + d0(i,j,k) * (s_now - s_prev) )/s_now;
          dmicro[i][j][k] = std::max(d_new, dmicro_min);
        }
      }

    } else {
      dmicro[i][j][k] = boil::unreal;
    }/* area vapor ? 0.0 */
  } /* ibody cells */

  dmicro.exchange();

  //real xxx;
  //std::cin>>xxx;

  return;
}
#endif
