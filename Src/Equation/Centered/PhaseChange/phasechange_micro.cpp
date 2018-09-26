#include "phasechange.h"
//#define ST_LENGTH
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange::micro(Vector * vec, const Scalar * diff_eddy) {

#ifdef DEBUG
  std::cout<<"pc.micro: "<<boil::cart.iam()<<"\n";
#endif

  if(nucl==NULL){
    boil::oout<<"# phasechange_micro:Error!  Nucleation is not defined.\n";
    exit(0);
  }

  /* initialize sum variables */
  real smdot_pos_macro_overwrite = 0.0;
  real smdot_neg_macro_overwrite = 0.0;
  for(int i=0; i<=6; i++){
    area_sum[i]=0.0;
    area_l[i]=0.0;
    area_v[i]=0.0;
    hflux_micro[i]=0.0;
    hflux_total[i]=0.0;
    hflux_vapor[i]=0.0;
  }

  /* initialize array */
#ifdef ST_MICRO
  txv =0.0;  // surface tension in x
  tyv =0.0;  // surface tension in y
  tzv =0.0;  // surface tension in z
#endif
  txl =0.0;  // normal vector x in micro layer
  tyl =0.0;  // normal vector y in micro layer
  tzl =0.0;  // normal vector z in micro layer

  /* store area of vapor to dSprev, if not stored */
  if (!nucl->store_dSprev) {
    str_dSprev();
  }

  /*------------------------------------+
  |  set nucleation site active or not  |
  +------------------------------------*/
  //nucl->st_active();

  /*---------------------+/
  |  intermediate dmicro  |
  +----------------------*/
  dmicro_intermediate();

  /*---------------------------------+
  |  calculate mdot on wall boundary |
  +---------------------------------*/
  for( int b=0; b<clr.bc().count(); b++ ) {
    if(clr.bc().type_decomp(b)) continue;
    if( clr.bc().type(b) == BndType::wall()) {
      int iof=0, jof=0, kof=0;
      Dir d = clr.bc().direction(b);
      int ndir = int(d);
      if(d != Dir::undefined()) {
        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;
        for_vijk( clr.bc().at(b), i,j,k ){
          // skip solid cells
          if(dom->ibody().off(i+iof,j+jof,k+kof)) continue;

          // only fluid cell comes here
          /*--------------------------------------------------+
          |  Note:  (i    , j    , k    ) is in solid domain  |
          |         (i+iof, j+jof, k+kof) is in fluid domain  |
          +--------------------------------------------------*/
          int ii=i+iof;
          int jj=j+jof;
          int kk=k+kof;
          real vol = dV(ii,jj,kk);

          /* area */
          real area = fabs(iof)*clr.dSx(ii,jj,kk)
                    + fabs(jof)*clr.dSy(ii,jj,kk)
                    + fabs(kof)*clr.dSz(ii,jj,kk);
          area_sum[ndir] += area;
          area_l[ndir] += area *clr[ii][jj][kk];
          area_v[ndir] += area *(1.0-clr[ii][jj][kk]);

          if (clr[ii][jj][kk] < 0.5) {
 
            /* micro region model */
            real dt = time->dt();
            if (nucl->dmicro[ii][jj][kk] >boil::mega) {
              nucl->dmicro[ii][jj][kk] = nucl->dmicro0(ii,jj,kk);
            }

            if (nucl->dmicro[ii][jj][kk] <= nucl->dmicro_min+boil::pico) {

              /* micro layer does not exist */
              /* heat flux for vapor cell */
              real lc=lambdav;
              if (diff_eddy) lc += (*diff_eddy)[ii][jj][kk]*cpv/rhov/turbP;
              real alen = fabs(tpr.xc(i)-tpr.xc(i+iof))
                        + fabs(tpr.yc(j)-tpr.yc(j+jof))
                        + fabs(tpr.zc(k)-tpr.zc(k+kof));
              hflux_total[ndir] += lc*(tpr[i][j][k]-tpr[i+iof][j+jof][k+kof])
                                /alen*area;
              hflux_vapor[ndir] += lc*(tpr[i][j][k]-tpr[i+iof][j+jof][k+kof])
                                /alen*area;

	    } else {

              /* micro layer exists */
	      real dmicro_new;
              //dmicro_new = nucl->dmicro[ii][jj][kk]
              //           - dt/rhol*lambdal/latent
              //           * ( tpr[i][j][k] - tsat ) / nucl->dmicro[ii][jj][kk];
              dmicro_new = nucl->dmicro[ii][jj][kk]
                         - dt/rhol/latent
                         * ( tpr[i][j][k] - tsat ) 
                         / (nucl->dmicro[ii][jj][kk]/lambdal + resint);
	      dmicro_new = max(dmicro_new, nucl->dmicro_min);

              // micro-layer thicnkess doesn't increase
              dmicro_new = min(dmicro_new,  nucl->dmicro[ii][jj][kk]); 

              phi[ii][jj][kk] = - rhol * (dmicro_new - nucl->dmicro[ii][jj][kk])
                              / dt * area / vol;
              clrs[ii][jj][kk] = -1.0/rhol * phi[ii][jj][kk];
	      nucl->dmicro[ii][jj][kk]=dmicro_new;

              /* heat flux */
              real Mdot = phi[ii][jj][kk] * vol / area;
              hflux_total[ndir] += Mdot * latent;
              hflux_micro[ndir] += Mdot * latent;

            }

          } else {

            /* heat flux for liquid cell */
            real lc=lambdal;
            if (diff_eddy) lc += (*diff_eddy)[ii][jj][kk]*cpl/rhol/turbP;
            real alen = fabs(tpr.xc(i)-tpr.xc(i+iof))
                      + fabs(tpr.yc(j)-tpr.yc(j+jof))
                      + fabs(tpr.zc(k)-tpr.zc(k+kof));
            hflux_total[ndir] += lc*(tpr[i][j][k]-tpr[i+iof][j+jof][k+kof])
                              /alen*area;

          }

          /*------------------+
          |  surface tension  |
          +------------------*/
          real nxx=nx[ii][jj][kk];
          real nyy=ny[ii][jj][kk];
          real nzz=nz[ii][jj][kk];
#ifdef ST_MICRO
          real drdx=(fluid()->rho(ii+1,jj,kk)-fluid()->rho(ii-1,jj,kk))
                   /(dxw(ii)+dxe(ii));
          real drdy=(fluid()->rho(ii,jj+1,kk)-fluid()->rho(ii,jj-1,kk))
                   /(dys(jj)+dyn(jj));
          real drdz=(fluid()->rho(ii,jj,kk+1)-fluid()->rho(ii,jj,kk-1))
                   /(dzb(kk)+dzt(kk));
          if(iof!=0){
            nxx=0.0;
            drdx=0.0;
          } else if(jof!=0){
            nyy=0.0;
            drdy=0.0;
          } else if(kof!=0){
            nzz=0.0;
            drdz=0.0;
          }
          normalize(nxx,nyy,nzz);
          real agradr=sqrt(drdx*drdx+drdy*drdy+drdz*drdz);
          real coefr =agradr*(fluid()->rho(ii,jj,kk))/(rhoave*rhodlt);
          txv[ii][jj][kk] = pcfrc(tpr[i][j][k]) * area * coefr * nxx;
          tyv[ii][jj][kk] = pcfrc(tpr[i][j][k]) * area * coefr * nyy;
          tzv[ii][jj][kk] = pcfrc(tpr[i][j][k]) * area * coefr * nzz;
#endif
          txl[ii][jj][kk] = -nxx;
          tyl[ii][jj][kk] = -nyy;
          tzl[ii][jj][kk] = -nzz;
        }
      }
    }
  }

#ifdef IB
  /*-------------------------------------+
  |  calculate mdot on immersed boundary |
  +-------------------------------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=dom->ibody().nwx(i,j,k);
    real uy=dom->ibody().nwy(i,j,k);
    real uz=dom->ibody().nwz(i,j,k);
    Dir d = Dir::undefined();
    if (fabs(uz)>0.707) {
      d = Dir::kmin();
      if(uz>0.707){
        d = Dir::kmax();
      }
    } else if (fabs(ux)>0.707) {
      d = Dir::imin();
      if(ux>0.707){
        d = Dir::imax();
      }
    } else if (fabs(uy)>0.707) {
      d = Dir::jmin();
      if(uy>0.707){
        d = Dir::jmax();
      }
    } else {
      std::cout<<"phasechange_micro: Underdevelopment!!!\n";
      exit(0);
    }

    int iof=0, jof=0, kof=0;
    int ndir = 6;  // ibody()
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*--------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in fluid domain  |
    |         (i+iof, j+jof, k+kof) is in solid domain  |
    +--------------------------------------------------*/
    real vol = dV(i,j,k);

    /* area */
    real area = fabs(iof)*clr.dSx(i,j,k)
              + fabs(jof)*clr.dSy(i,j,k)
              + fabs(kof)*clr.dSz(i,j,k);
    real a_vapor = nucl->area_vapor(i,j,k,d);
    real clrc = min(1.0,max(0.0,clr[i][j][k]));
    bool inclNoInterface = incl_no_interface(i,j,k);
    if (incl_no_interface(i,j,k)) {
      if (clrc>0.5) {
        clrc=1.0;
      } else {
        clrc=0.0;
      }
    }
    area_sum[ndir] += area;
    area_l[ndir] += area - a_vapor;
    area_v[ndir] += a_vapor;

    int ii=i+iof;
    int jj=j+jof;
    int kk=k+kof;
    real alen = fabs(tpr.xc(i)-tpr.xc(ii))
              + fabs(tpr.yc(j)-tpr.yc(jj))
              + fabs(tpr.zc(k)-tpr.zc(kk));

    real dw = 0.5 * tpr.dzc(kk);  // half cell size in wall
    real df = 0.5 * tpr.dzc(k);   // half cell size in fluid
    real lambdas = solid()->lambda(ii,jj,kk);  // lambda solid

    /* cell includes vapor: includes micro-layer */
    if (a_vapor == 0.0) {

      /* liquid cell */
      nucl->dmicro[i][j][k]=boil::exa;

      /* heat flux for liquid cell */
      real laml = lambdal;
      if (diff_eddy) laml += (*diff_eddy)[i][j][k]*cpl/rhol/turbP; 
      real tw = ( df * lambdas * tpr[ii][jj][kk]
                + dw * laml    * tpr[i ][j ][k ] )
                / ( df * lambdas + dw * laml);
      hflux_total[ndir] += lambdas*(tpr[ii][jj][kk]-tw)/dw*area;

    } else if ( approx (a_vapor, area, area*boil::micro)) {

      /* vapor cell */

      if ( nucl->dmicro[i][j][k] <= nucl->dmicro_min+boil::pico // depleted
        || nucl->dmicro[i][j][k] >  boil::mega){                // full vapor

        /* micro layer depleted or full vapor cell */
        /* temperature at heat transfer surface */
        real lamv = lambdav;
        if (diff_eddy) lamv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
        real tw = ( df * lambdas * tpr[ii][jj][kk]
                  + dw * lamv    * tpr[i ][j ][k ] )
                  / ( df * lambdas + dw * lamv);

        /* heat flux for vapor cell */
        hflux_total[ndir] += lambdas*(tpr[ii][jj][kk]-tw)/dw*area;
        hflux_vapor[ndir] += lambdas*(tpr[ii][jj][kk]-tw)/dw*area;

      } else {

        /* micro layer exists */
        real dmicro_new;
        real dt = time->dt();

        /* temperature at heat transfer surface */
        //real tw = ( nucl->dmicro[i][j][k] * lambdas * tpr[ii][jj][kk]
        //          + dw * lambdal * tsat )
        //          / ( nucl->dmicro[i][j][k] * lambdas + dw * lambdal);
        real qtmp = (tpr[ii][jj][kk] - tsat)
                  / ( dw/lambdas + nucl->dmicro[i][j][k]/lambdal + resint); 

        dmicro_new = nucl->dmicro[i][j][k] - dt / rhol * qtmp / latent;
        dmicro_new = max(dmicro_new, nucl->dmicro_min);

        // micro-layer thicnkess doesn't increase
        dmicro_new = min(dmicro_new,  nucl->dmicro[i][j][k]); 

        if (phi[i][j][k]>=0.0) {
          smdot_pos_macro_overwrite += phi[i][j][k]*vol;
        } else {
          smdot_neg_macro_overwrite += phi[i][j][k]*vol;
        }
        /* overwrite phi */
        phi[i][j][k] = - rhol * (dmicro_new - nucl->dmicro[i][j][k])
                       / dt * a_vapor / vol;
        clrs[i][j][k] = -1.0/rhol * phi[i][j][k];
        tprs[i+iof][j+jof][k+kof]  -= vol*phi[i][j][k]*latent;
        tprs[i][j][k] -= cpv *(tpr[i][j][k]-tsat)*(1.0/rhov-1.0/rhol)
                         *phi[i][j][k]*vol;
        nucl->dmicro[i][j][k]=dmicro_new;

        /* heat flux */
        // it doesn't take into account max(dmicro_new, nucl->dmicro_min)
        hflux_total[ndir] += qtmp*a_vapor;
        hflux_micro[ndir] += qtmp*a_vapor;
        // it does take into account max(dmicro_new, nucl->dmicro_min)
        //hflux_total[ndir] += phi[i][j][k] * vol * latent;
        //hflux_micro[ndir] += phi[i][j][k] * vol * latent;

      }

    } else {

      /* liquid and vapor cell */

      /* temperature at heat transfer surface covered with liquid */
      real laml = lambdal;
      if (diff_eddy) laml += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      real twl = ( df * lambdas * tpr[ii][jj][kk]
                 + dw * laml    * tpr[i ][j ][k ] )
                 / ( df * lambdas + dw * laml);

      if ( nucl->dmicro[i][j][k] <= nucl->dmicro_min+boil::pico // depleted
        || nucl->dmicro[i][j][k] >  boil::mega){                // full vapor

        /* micro layer depleted or full vapor cell */
        /* temperature at heat transfer surface covered with vapor */
        real lamv = lambdav;
        if (diff_eddy) lamv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
        real twv = ( df * lambdas * tpr[ii][jj][kk] 
                   + dw * lamv    * tpr[i ][j ][k ] )
                   / ( df * lambdas + dw * lamv);
        /* heat flux for vapor cell */
        hflux_total[ndir] += lambdas*(tpr[ii][jj][kk]-twv)/dw*a_vapor
                           + lambdas*(tpr[ii][jj][kk]-twl)/dw*(area-a_vapor);
        hflux_vapor[ndir] += lambdas*(tpr[ii][jj][kk]-twv)/dw*a_vapor;

      } else {

        /* micro layer exists */
        real dmicro_new;
        real dt = time->dt();

        /* temperature at heat transfer surface covered with micro-layer */
        //real twv = ( nucl->dmicro[i][j][k] * lambdas * tpr[ii][jj][kk]
        //           + dw * lambdal * tsat )
        //           / ( nucl->dmicro[i][j][k] * lambdas + dw * lambdal);

        real qtmp = (tpr[ii][jj][kk] - tsat)
                  / ( dw/lambdas + nucl->dmicro[i][j][k]/lambdal + resint);
        
        dmicro_new = nucl->dmicro[i][j][k] - dt / rhol * qtmp / latent;
        dmicro_new = max(dmicro_new, nucl->dmicro_min);
       
        // micro-layer thicnkess doesn't increase
        dmicro_new = min(dmicro_new,  nucl->dmicro[i][j][k]); 
 
        if (phi[i][j][k]>=0.0) {
          smdot_pos_macro_overwrite += phi[i][j][k]*vol;
        } else {
          smdot_neg_macro_overwrite += phi[i][j][k]*vol;
        }
        /* overwrite phi */
        phi[i][j][k] = - rhol * (dmicro_new - nucl->dmicro[i][j][k])
                       / dt * a_vapor / vol;
        clrs[i][j][k] = -1.0/rhol * phi[i][j][k];
        tprs[i+iof][j+jof][k+kof]  -= vol*phi[i][j][k]*latent;
        tprs[i][j][k] -= cpv *(tpr[i][j][k]-tsat)*(1.0/rhov-1.0/rhol)
                         *phi[i][j][k]*vol;
        nucl->dmicro[i][j][k]=dmicro_new;

        /* heat flux */
        // it doesn't take into account max(dmicro_new, nucl->dmicro_min)
        hflux_total[ndir] += qtmp*a_vapor
                             + lambdas*(tpr[ii][jj][kk]-twl)/dw*(area-a_vapor);
        hflux_micro[ndir] += qtmp*a_vapor;
        // it does take into account max(dmicro_new, nucl->dmicro_min)
        //hflux_total[ndir] += phi[i][j][k] * vol * latent
        //                   + lambdas*(tpr[ii][jj][kk]-twl)/dw*(area-a_vapor);
        //hflux_micro[ndir] += phi[i][j][k] * vol * latent;

      }

    }

      /* surface tension */ 
      real nxx=nx[i][j][k];
      real nyy=ny[i][j][k];
      real nzz=nz[i][j][k];
      real drdx=(fluid()->rho(i+1,j,k)-fluid()->rho(i-1,j,k))
               /(dxw(i)+dxe(i));
      real drdy=(fluid()->rho(i,j+1,k)-fluid()->rho(i,j-1,k))
               /(dys(j)+dyn(j));
      real drdz=(fluid()->rho(i,j,k+1)-fluid()->rho(i,j,k-1))
               /(dzb(k)+dzt(k));
      if(iof!=0){
        nxx=0.0;
        drdx=0.0;
      } else if(jof!=0){
        nyy=0.0;
        drdy=0.0;
      } else if(kof!=0){
        nzz=0.0;
        drdz=0.0;
      }
      normalize(nxx,nyy,nzz);
#ifdef ST_MICRO
      real agradr=sqrt(drdx*drdx+drdy*drdy+drdz*drdz);
      real coefr =agradr*(fluid()->rho(i,j,k))/(rhoave*rhodlt);
      real tprw = tpr[i+iof][j+jof][k+kof];
      txv[i][j][k] = pcfrc(tprw) * area * coefr * nxx;
      tyv[i][j][k] = pcfrc(tprw) * area * coefr * nyy;
      tzv[i][j][k] = pcfrc(tprw) * area * coefr * nzz;
#endif
      txl[i][j][k] = -nxx;
      tyl[i][j][k] = -nyy;
      tzl[i][j][k] = -nzz;
  }
#endif

  phi.exchange_all();
#ifdef ST_MICRO
  txv.exchange_all();
  tyv.exchange_all();
  tzv.exchange_all();
#endif
  txl.exchange();
  tyl.exchange();
  tzl.exchange();

  boil::cart.sum_real_n(area_sum,7);
  boil::cart.sum_real_n(area_l,7);
  boil::cart.sum_real_n(area_v,7);
  boil::cart.sum_real_n(hflux_total,7);
  boil::cart.sum_real_n(hflux_micro,7);
  boil::cart.sum_real_n(hflux_vapor,7);

  /*----------------------------------------+
  |  surface tension of micro region model  |
  +----------------------------------------*/
#ifdef ST_MICRO
  Comp m;
  m = Comp::u();
  buff=0.0;
  for_ijk(i,j,k) {
    real coef=1.0;
    //if(stmp[i][j][k]==0.0) coef=0.0;
    (*vec)[m][i  ][j][k] += 0.5 * txv[i][j][k] * coef;
    (*vec)[m][i+1][j][k] += 0.5 * txv[i][j][k] * coef;

    if(i==si()) buff[si()][j][k] = 0.5 * txv[i][j][k] * coef;
    if(i==ei()) buff[ei()][j][k] = 0.5 * txv[i][j][k] * coef;
  }

  buff.exchange(0);
  for_vmijk((*vec),m,i,j,k) {
    if(i==si())   (*vec)[m][i][j][k] += buff[si()-1][j][k];
    if(i==ei()+1) (*vec)[m][i][j][k] += buff[ei()+1][j][k];
  }

  m = Comp::v();
  buff=0.0;
  for_ijk(i,j,k) {
    real coef=1.0;
    //if(stmp[i][j][k]==0.0) coef=0.0;
    (*vec)[m][i][j  ][k] += 0.5 * tyv[i][j][k] * coef;
    (*vec)[m][i][j+1][k] += 0.5 * tyv[i][j][k] * coef;

    if(j==sj()) buff[i][sj()][k] = 0.5 * tyv[i][j][k] * coef;
    if(j==ej()) buff[i][ej()][k] = 0.5 * tyv[i][j][k] * coef;
  }

  buff.exchange(1);
  for_vmijk((*vec),m,i,j,k) {
    if(j==sj())   (*vec)[m][i][j][k] += buff[i][sj()-1][k];
    if(j==ej()+1) (*vec)[m][i][j][k] += buff[i][ej()+1][k];
  }

  m = Comp::w();
  buff=0.0;
  for_ijk(i,j,k) {
    real coef=1.0;
    //if(stmp[i][j][k]==0.0) coef=0.0;
    (*vec)[m][i][j][k  ] += 0.5 * tzv[i][j][k] * coef;
    (*vec)[m][i][j][k+1] += 0.5 * tzv[i][j][k] * coef;

    if(k==sk()) buff[i][j][sk()] = 0.5 * tzv[i][j][k] * coef;
    if(k==ek()) buff[i][j][ek()] = 0.5 * tzv[i][j][k] * coef;
  }
  buff.exchange(2);
  for_vmijk((*vec),m,i,j,k) {
    if(k==sk())   (*vec)[m][i][j][k] += buff[i][j][sk()-1];
    if(k==ek()+1) (*vec)[m][i][j][k] += buff[i][j][ek()+1];
  }
#endif

  /* shift mdot */
  //micro_shift();

#ifdef DEBUG
  boil::plot->plot(clr, tpr, phi, "clr-tpr-mdot",  time->current_step());
#endif

  /*-------------------------+
  |  calculate source terms  |
  +-------------------------*/
  /* clrs is already computed in this function */
  sources_fext();
  //sources_tprs();
  sources_sum();

  boil::cart.sum_real(&smdot_pos_macro_overwrite);
  boil::cart.sum_real(&smdot_neg_macro_overwrite);

  boil::oout<<"phasechange_micro: time= "<<time->current_time()
            <<" smdot_pos_macro= "<<smdot_pos_macro-smdot_pos_macro_overwrite
            <<" smdot_pos_incl.micro= "<<smdot_pos
            <<" smdot_neg_macro= "<<smdot_neg_macro-smdot_neg_macro_overwrite
            <<"\n";

#ifdef DEBUG
  boil::plot->plot(clr, nucl->dmicro, phi, clrs, 
                  "clr-dmicro-mdot-clrs",  time->current_step());
  exit(0);
#endif

  /* store area of vapor to dSprev */
  str_dSprev();

}

/*-----------------------------------------------------------------------------+
 '$Id: phasechange_micro.cpp,v 1.19 2016/04/15 09:27:50 sato Exp $'/
+-----------------------------------------------------------------------------*/
