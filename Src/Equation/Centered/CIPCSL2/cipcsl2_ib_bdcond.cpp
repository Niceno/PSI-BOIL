#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::ib_bdcond(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for immersed boundary.
******************************************************************************/

  if(sca.domain()->ibody().ncall() == 0) return;

  /* cell to face */
#if 0
  for(int cc=0; cc<dom->ibody().nccells(); cc++){

    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    /* set direction */    // crude code!!!
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
      std::cout<<"cipcsl2_ib_insert_bc: Underdevelopment!!!\n";
      exit(0);
    }

    real area;
    Comp m;
    if(d == Dir::imin() || d == Dir::imax()) { 
      m = Comp::u();
      area = sca.dyc(j)*sca.dzc(k);
    } else if(d == Dir::jmin() || d == Dir::jmax()) { 
      m = Comp::v();
      area = sca.dxc(i)*sca.dzc(k);
    } else { 
      m = Comp::w();
      area = sca.dxc(i)*sca.dyc(j);
    }

    int iof=0, jof=0, kof=0;
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*----------------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in fluid domain  |
    |         (i+iof, j+jof, k+kof) is in wall                  |
    +----------------------------------------------------------*/
    if(dom->ibody().off(i+iof,j+jof,k+kof)){
      iof=max(0,iof);
      jof=max(0,jof);
      kof=max(0,kof);
      sxyz[m][i+iof][j+jof][k+kof] = sca[i][j][k]*area;
    }
  }
  sxyz.exchange_all();

  /* face to edge */

  for(int cc=0; cc<dom->ibody().nccells(); cc++){

    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    /* set direction */    // crude code!!!
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
      std::cout<<"cipcsl2_ib_insert_bc: Underdevelopment!!!\n";
      exit(0);
    }

    real area;
    Comp m;
    if(d == Dir::imin() || d == Dir::imax()) { 
      m = Comp::u();
    } else if(d == Dir::jmin() || d == Dir::jmax()) { 
      m = Comp::v();
    } else { 
      m = Comp::w();
    }

    int iof=0, jof=0, kof=0;
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*----------------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in fluid domain  |
    |         (i+iof, j+jof, k+kof) is in wall                  |
    +----------------------------------------------------------*/
    if(dom->ibody().off(i+iof,j+jof,k+kof)){
      if(d != Dir::kmin()){
        boil::oout<<"cipcsl2_ib_bdcond: Not implemented yet !!!\n";
        boil::oout<<"Dir= "<<d<<"\n";
        exit(0);
      }
      if(d == Dir::kmin()){
        if(sca.dyc(j-1)==0.0) {
          scheme.sigx[i][j  ][k] = 
            ( sxyz[m][i][j  ][k]/(sca.dxc(i)*sca.dyc(j))) * sca.dxc(i);
        } else {
          scheme.sigx[i][j  ][k] = 
            ( sxyz[m][i][j-1][k]/(sca.dxc(i)*sca.dyc(j-1))
            + sxyz[m][i][j  ][k]/(sca.dxc(i)*sca.dyc(j)) )
            * 0.5 * sca.dxc(i);
        }
        if(sca.dyc(j+1)==0.0) {
          scheme.sigx[i][j+1][k] = 
            ( sxyz[m][i][j  ][k]/(sca.dxc(i)*sca.dyc(j)) ) * sca.dxc(i);
        } else {
          scheme.sigx[i][j+1][k] = 
            ( sxyz[m][i][j+1][k]/(sca.dxc(i)*sca.dyc(j+1))
            + sxyz[m][i][j  ][k]/(sca.dxc(i)*sca.dyc(j)) )
            * 0.5 * sca.dxc(i);
        }
        if(sca.dxc(i-1)==0) {
          scheme.sigy[i  ][j][k] = 
            ( sxyz[m][i  ][j][k]/(sca.dxc(i  )*sca.dyc(j))) * sca.dyc(j);
        } else {
          scheme.sigy[i  ][j][k] = 
            ( sxyz[m][i-1][j][k]/(sca.dxc(i-1)*sca.dyc(j))
            + sxyz[m][i  ][j][k]/(sca.dxc(i  )*sca.dyc(j)))
            * 0.5 * sca.dyc(j);
        }
        if(sca.dxc(i+1)==0) {
          scheme.sigy[i+1][j][k] = 
            ( sxyz[m][i  ][j][k]/(sca.dxc(i  )*sca.dyc(j))) * sca.dyc(j);
        } else {
          scheme.sigy[i+1][j][k] = 
            ( sxyz[m][i+1][j][k]/(sca.dxc(i+1)*sca.dyc(j))
            + sxyz[m][i  ][j][k]/(sca.dxc(i  )*sca.dyc(j)))
            * 0.5 * sca.dyc(j);
        }
        for(int ii=0; ii<=1; ii++){
          for(int jj=0; jj<=1; jj++){
            int it=i+ii;
            int jt=j+jj;
            int icount = 0;
            real app = sca.dxc(it  )*sca.dyc(jt  );
            real amp = sca.dxc(it-1)*sca.dyc(jt  );
            real apm = sca.dxc(it  )*sca.dyc(jt-1);
            real amm = sca.dxc(it-1)*sca.dyc(jt-1);
            scheme.f[it][jt][k] = 0.0;
            if(app!=0.0) {
              scheme.f[it][jt][k] += sxyz[m][it  ][jt  ][k]/app;
              icount ++;
            }
            if(amp!=0.0) {
              scheme.f[it][jt][k] += sxyz[m][it-1][jt  ][k]/amp;
              icount ++;
            }
            if(apm!=0.0) {
              scheme.f[it][jt][k] += sxyz[m][it  ][jt-1][k]/apm;
              icount ++;
            }
            if(amm!=0.0) {
              scheme.f[it][jt][k] += sxyz[m][it-1][jt-1][k]/amm;
              icount ++;
            }
	    if(icount==0){
		    cout<<"Error!!!\n";
		    exit(0);
            } else {
              scheme.f[it][jt][k] /= real(icount);
            }
            //scheme.f[it][jt][k] =
            // ( sxyz[m][it  ][jt  ][k]/(sca.dxc(it  )*sca.dyc(jt  ))
            // + sxyz[m][it-1][jt  ][k]/(sca.dxc(it-1)*sca.dyc(jt  ))
            // + sxyz[m][it  ][jt-1][k]/(sca.dxc(it  )*sca.dyc(jt-1))
            // + sxyz[m][it-1][jt-1][k]/(sca.dxc(it-1)*sca.dyc(jt-1)))*0.25;
          }
        }
      }
    }
  }

  scheme.f.exchange_all();
  scheme.sigx.exchange_all();
  scheme.sigy.exchange_all();
  scheme.sigz.exchange_all();
#endif
}
