#include "model.h"

/***************************************************************************//**
*  Computes eddy viscosity using Smagorinsky SGS model. Unit is [kg/(ms)]. 
*******************************************************************************/
void Model::smagorinsky( const Momentum * mom, Scalar * mu_t, real c_s ) const {

  /* take handy aliases */
  const Vector & uvw   =   mom->val();
  const Matter & fluid = * mom->fluid();

  assert(mu_t->domain() == uvw.domain() );

  for_vijk( (*mu_t), i, j, k) {
    /*----------------------------------------------+
    |  shear = sqrt( 2 * s_ij * s_ij )       [1/s]  |
    |  s_ij = 1/2 ( du_i/dx_j + du_j/dx_i )  [1/s]  |
    +----------------------------------------------*/
    real shear = uvw.shear_magn(i,j,k);

    /*------------------------------------+
    |  delta = (dx * dy * dz)^(1/3)  [m]  |
    +------------------------------------*/
    real delta = pow((mu_t->dxc(i)*mu_t->dyc(j)*mu_t->dzc(k)),1.0/3.0); 

    /*-----------------------------------------------+
    |  mu_t = rho * (Cs * delta)^2 * |s|  [kg/(ms)]  |
    +-----------------------------------------------*/
    (*mu_t)[i][j][k] = fluid.rho(i,j,k) * pow((c_s*delta),2.0)*shear; 
  }

  /*---------------------------------------------------+
  |  eliminate excessive viscosity close to the wall.  |
  |  it will be handled by wall function in any case.  |
  +---------------------------------------------------*/
  if( uvw.bc(Comp::u()).type_here( Dir::imin(), BndType::wall()) ) {
    const int i=mu_t->si();
    for_vjk( (*mu_t),j,k ) (*mu_t)[i][j][k] = (*mu_t)[i+1][j][k];
  }

  if( uvw.bc(Comp::u()).type_here( Dir::imax(), BndType::wall()) ) {
    const int i=mu_t->ei();
    for_vjk( (*mu_t),j,k ) (*mu_t)[i][j][k] = (*mu_t)[i-1][j][k];
  }

  if( uvw.bc(Comp::u()).type_here( Dir::jmin(), BndType::wall()) ) {
    const int j=mu_t->sj();
    for_vik( (*mu_t),i,k ) (*mu_t)[i][j][k] = (*mu_t)[i][j+1][k];
  }

  if( uvw.bc(Comp::u()).type_here( Dir::jmax(), BndType::wall()) ) {
    const int j=mu_t->ej();
    for_vik( (*mu_t),i,k ) (*mu_t)[i][j][k] = (*mu_t)[i][j-1][k];
  }

  if( uvw.bc(Comp::u()).type_here( Dir::kmin(), BndType::wall()) ) {
    const int k=mu_t->sk();
    for_vij( (*mu_t),i,j ) (*mu_t)[i][j][k] = (*mu_t)[i][j][k+1];
  }

  if( uvw.bc(Comp::u()).type_here( Dir::kmax(), BndType::wall()) ) {
    const int k=mu_t->ek();
    for_vij( (*mu_t),i,j ) (*mu_t)[i][j][k] = (*mu_t)[i][j][k-1];
  }

  /* for immersed boundary */
#if 0
  /* take handy aliases */
  const Body   & ibody =   mom->domain()->ibody(); 
  if( ibody.ncall() > 0 ) {
    for(int cc=0; cc<ibody.nccells(); cc++){
      int i,j,k;
      ibody.ijk(cc,&i,&j,&k);

      real sum   = 0.0;
      int  nwall = 0;

      if( ibody.off(i-1,j,k) ) {
        sum += (*mu_t)[i+1][j][k];
        nwall++;
      }
      if( ibody.off(i+1,j,k) ) {
        sum += (*mu_t)[i-1][j][k];
        nwall++;
      }

      if( ibody.off(i,j-1,k) ) {
        sum += (*mu_t)[i][j+1][k];
        nwall++;
      }

      if( ibody.off(i,j+1,k) ) {
        sum += (*mu_t)[i][j-1][k];
        nwall++;
      }

      if( ibody.off(i,j,k-1) ) {
        sum += (*mu_t)[i][j][k+1];
        nwall++;
      }

      if( ibody.off(i,j,k+1) ) {
        sum += (*mu_t)[i][j][k-1];
        nwall++;
      }

      if (nwall>=1) {
        (*mu_t)[i][j][k] = sum/real(nwall);
      } else {
        std::cout<<"smagorinsky: nwall=0 at proc= "<<boil::cart.iam()<<" i= "
                 <<i<<" j= "<<j<<" k= "<<k<<"\n";
      }
    }
  }
#endif

  /*----------------------------------+ 
  |  exchange buffers. important to   |
  |  exchange in the corners as well  |
  +----------------------------------*/
  mu_t->exchange_all();
}

/*-----------------------------------------------------------------------------+
 '$Id: model_smagorinsky.cpp,v 1.9 2015/08/10 12:30:25 sato Exp $'/
+-----------------------------------------------------------------------------*/
