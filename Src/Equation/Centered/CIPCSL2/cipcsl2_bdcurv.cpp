#include "cipcsl2.h"
#include <iomanip>
using namespace std;

real divnorm(const Scalar & nx, const Scalar & ny,const Scalar & nz
            ,const int & i, const int & j, const int & k);

/******************************************************************************/
void CIPCSL2::bdcurv(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate curvature on wall boundary condition.
*         Reference: J.U.Brackbill, et al.,"A Continum method for modeling
*                    surface tension",J.Comp.phys.,Vol.100,pp.335-354,1992
*                    Equation (41)
*******************************************************************************/
  boil::timer.start("cipcsl2 bdcurv");

  /*----------------------------------------+
  |  calculate normal vector on node point  |
  +----------------------------------------*/

  /* cal normal vector on node point */
  gradphi(sca);

  /* cal normal vector on wall considering contact angle */
  wall_norm(sca);


#ifdef IB
  /* cal normal vector on immersed boundary considering contact angle */
  if(dom->ibody().nccells() > 0) {
    ib_norm(sca);
  }
#endif

  /*---------------------------------------------+
  |  calculate curvature. curvature=-div(norm)   |
  +---------------------------------------------*/
  if(sca.bc().type_here(Dir::imin(), BndType::wall())) {
    int ii=si();
    for_vjk(kappa,j,k) {
      kappa[ii][j][k]=-divnorm(nx,ny,nz,ii,j,k);
      kappa[ii-1][j][k]=kappa[ii][j][k];
    }
  }

  if(sca.bc().type_here(Dir::imax(), BndType::wall())) {
    int ii=ei();
    for_vjk(kappa,j,k) {
      kappa[ii][j][k]=-divnorm(nx,ny,nz,ii,j,k);
      kappa[ii+1][j][k]=kappa[ii][j][k];
    }
  }

  if(sca.bc().type_here(Dir::jmin(), BndType::wall())) {
    int jj=sj();
    for_vik(kappa,i,k) {
      kappa[i][jj][k]=-divnorm(nx,ny,nz,i,jj,k);
      kappa[i][jj-1][k]=kappa[i][jj][k];
    }
  }

  if(sca.bc().type_here(Dir::jmax(), BndType::wall())) {
    int jj=ej();
    for_vik(kappa,i,k) {
      kappa[i][jj][k]=-divnorm(nx,ny,nz,i,jj,k);
      kappa[i][jj+1][k]=kappa[i][jj][k];
    }
  }

  if(sca.bc().type_here(Dir::kmin(), BndType::wall())) {
    int kk=sk();
    for_vij(kappa,i,j) {
      kappa[i][j][kk]=-divnorm(nx,ny,nz,i,j,kk);
      kappa[i][j][kk-1]=kappa[i][j][kk];
    }
  }

  if(sca.bc().type_here(Dir::kmax(), BndType::wall())) {
    int kk=ek();
    for_vij(kappa,i,j) {
      kappa[i][j][kk]=-divnorm(nx,ny,nz,i,j,kk);
      kappa[i][j][kk+1]=kappa[i][j][kk];
    }
  }

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    kappa[i][j][k] = -divnorm(nx,ny,nz,i,j,k);
  }
#endif

  insert_bc_kappa(kappa);
  kappa.exchange();

#if 1
  /*------------------------------------+
  |  set_wflag for wall adjacent cells  |
  +------------------------------------*/
  set_wflag2();

  /*-------------------------------------------+
  |  interpolate curvature in interface cells  |
  +-------------------------------------------*/
  bdcurv_interface();

#if 0
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa_afterbdInter",
                     time->current_step());
#endif

  /*------------------------------------------------------+
  |  extrapolate curvature from interface cell to others  |
  |   in wall adjacent layer                              |
  +------------------------------------------------------*/
  bdcurv_interface_ext();


#if 0
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa_afterbdInterExt",
                     time->current_step());
#endif
#endif

  boil::timer.stop("cipcsl2 bdcurv");
  return;
}

/******************************************************************************/
real divnorm(const Scalar & nx, const Scalar & ny,const Scalar & nz
            ,const int & i, const int & j, const int & k){

  real nw = 0.25* (nx[i  ][j][k  ]+nx[i  ][j+1][k  ]
                  +nx[i  ][j][k+1]+nx[i  ][j+1][k+1]);
  real ne = 0.25* (nx[i+1][j][k  ]+nx[i+1][j+1][k  ]
                  +nx[i+1][j][k+1]+nx[i+1][j+1][k+1]);
  real ns = 0.25* (ny[i][j  ][k  ]+ny[i+1][j  ][k  ]
                  +ny[i][j  ][k+1]+ny[i+1][j  ][k+1]);
  real nn = 0.25* (ny[i][j+1][k  ]+ny[i+1][j+1][k  ]
                  +ny[i][j+1][k+1]+ny[i+1][j+1][k+1]);
  real nb = 0.25* (nz[i][j  ][k  ]+nz[i+1][j  ][k  ]
                  +nz[i][j+1][k  ]+nz[i+1][j+1][k  ]);
  real nt = 0.25* (nz[i][j  ][k+1]+nz[i+1][j  ][k+1]
                  +nz[i][j+1][k+1]+nz[i+1][j+1][k+1]);

  return( ((ne-nw)/(nx.dxc(i))
          +(nn-ns)/(ny.dyc(j))
          +(nt-nb)/(nz.dzc(k))) );

}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_bdcurv.cpp,v 1.5 2015/05/05 15:18:24 sato Exp $'/
+-----------------------------------------------------------------------------*/

