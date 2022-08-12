#include "vof.h"
#include <iomanip>

/******************************************************************************/
void VOF::ib_norm(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector on immersed boundary
*******************************************************************************/

  /*--------------------------------+
  |  normal vector of free surface  |
  +--------------------------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    if(dom->ibody().fPmmm(i,j,k)==0) ib_norm_cal(cc, i  ,j  ,k  );
    if(dom->ibody().fPpmm(i,j,k)==0) ib_norm_cal(cc, i+1,j  ,k  );
    if(dom->ibody().fPmpm(i,j,k)==0) ib_norm_cal(cc, i  ,j+1,k  );
    if(dom->ibody().fPppm(i,j,k)==0) ib_norm_cal(cc, i+1,j+1,k  );
    if(dom->ibody().fPmmp(i,j,k)==0) ib_norm_cal(cc, i  ,j  ,k+1);
    if(dom->ibody().fPpmp(i,j,k)==0) ib_norm_cal(cc, i+1,j  ,k+1);
    if(dom->ibody().fPmpp(i,j,k)==0) ib_norm_cal(cc, i  ,j+1,k+1);
    if(dom->ibody().fPppp(i,j,k)==0) ib_norm_cal(cc, i+1,j+1,k+1);
  }

#if 0
  boil::plot->plot(sca,nx,ny,nz, 
              "ib_norm-dist-nx-ny-nz", time->current_step());
#endif

  return;
}

/******************************************************************************/
void VOF::ib_norm_cal(const int cc, const int i, const int j, const int k){

    /*--------------------------------+
    |  normal vector of free surface  |
    +--------------------------------*/
    real nfx, nfy, nfz;
    nfx = nx[i][j][k];
    nfy = ny[i][j][k];
    nfz = nz[i][j][k];

    /*------------------------+
    |  normal vector of wall  |
    +------------------------*/
    real nwx, nwy, nwz;
    nwx =-0.25*((dom->ibody().dist(i,j  ,k  ) - dom->ibody().dist(i-1,j  ,k  ))
               +(dom->ibody().dist(i,j-1,k  ) - dom->ibody().dist(i-1,j-1,k  ))
               +(dom->ibody().dist(i,j  ,k-1) - dom->ibody().dist(i-1,j  ,k-1))
               +(dom->ibody().dist(i,j-1,k-1) - dom->ibody().dist(i-1,j-1,k-1)))
        / dxw(i);
    nwy =-0.25*((dom->ibody().dist(i-1,j,k  ) - dom->ibody().dist(i-1,j-1,k  ))
               +(dom->ibody().dist(i  ,j,k  ) - dom->ibody().dist(i  ,j-1,k  ))
               +(dom->ibody().dist(i-1,j,k-1) - dom->ibody().dist(i-1,j-1,k-1))
               +(dom->ibody().dist(i  ,j,k-1) - dom->ibody().dist(i  ,j-1,k-1)))
        / dys(j);
    nwz =-0.25*((dom->ibody().dist(i-1,j-1,k) - dom->ibody().dist(i-1,j-1,k-1))
               +(dom->ibody().dist(i  ,j-1,k) - dom->ibody().dist(i  ,j-1,k-1))
               +(dom->ibody().dist(i-1,j  ,k) - dom->ibody().dist(i-1,j  ,k-1))
               +(dom->ibody().dist(i  ,j  ,k) - dom->ibody().dist(i  ,j  ,k-1)))
        / dzb(k);
    normalize(nwx, nwy, nwz);
    //std::cout<<nwx<<" "<<nwy<<" "<<nwz<<"\n";

    real nout[3];
    nib(nwx,nwy,nwz,nfx,nfy,nfz,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];

  return;
}
/******************************************************************************/
void VOF::nib(const real & nwx, const real & nwy, const real & nwz
       , const real & nfx, const real & nfy, const real & nfz
       , real nout[]){
  // normal vector to wall:         (nwx, nwy, nwz)
  // normal vector to free surface: (nfx, nfy, nfz)
 
  real pin, ntanx, ntany, ntanz;

  // inner product of (nw.nf)
  pin = nwx*nfx + nwy*nfy + nwz*nfz;

  // nt_vector = nf_vector - (nf.nf)*nw_vector
  ntanx = nfx - pin * nwx;
  ntany = nfy - pin * nwy;
  ntanz = nfz - pin * nwz;
  normalize(ntanx, ntany, ntanz);
   
  // n_vector = nw_vector*cos(c) + nt_vector*sin(c): Eq.53
  real nnx = nwx*cos(cangle0) + ntanx*sin(cangle0);
  real nny = nwy*cos(cangle0) + ntany*sin(cangle0);
  real nnz = nwz*cos(cangle0) + ntanz*sin(cangle0);
  normalize(nnx,nny,nnz);
  nout[0] = nnx;
  nout[1] = nny;
  nout[2] = nnz;

}

