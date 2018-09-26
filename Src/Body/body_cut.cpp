#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::cut(const Scalar & phi) {

  phi.exchange_all();

  /*------------------------------------------------------------+ 
  |  define and compute scalar for storing nodal values of phi  |
  +------------------------------------------------------------*/
  Scalar nodal_phi( * phi.domain());

  for(int i=1; i<phi.ni(); i++)
    for(int j=1; j<phi.nj(); j++)
      for(int k=1; k<phi.nk(); k++) {
        nodal_phi[i][j][k] = .125 * (phi[i][j]  [k]   + phi[i-1][j]  [k]   +
                                   + phi[i][j-1][k]   + phi[i-1][j-1][k]   +
                                   + phi[i][j]  [k-1] + phi[i-1][j]  [k-1] +
                                   + phi[i][j-1][k-1] + phi[i-1][j-1][k-1]);  
      }

  /*--------------------------------------------------------+ 
  |  new polygon coordinates (size 6 should be enough, but  |
  |  in case cell lies on the body one might get 8 nodes)   |
  +--------------------------------------------------------*/
  real xp[8];
  real yp[8];
  real zp[8];

  /*---------------------------+ 
  |                            |
  |  browse through all cells  |
  |                            |
  +---------------------------*/
  for_vijk(phi,i,j,k) {

    /*-----------------------------------------------------------+
    |  compute normal as derivative of nodal indicator function  |
    +-----------------------------------------------------------*/
    real nx = 0.0;
    for(int jj=0; jj<2; jj++) 
      for(int kk=0; kk<2; kk++) 
        nx += (nodal_phi[i+1][j+jj][k+kk] - nodal_phi[i][j+jj][k+kk]); 
    nx *= (0.25 / phi.dxc(i));

    real ny = 0.0;
    for(int ii=0; ii<2; ii++) 
      for(int kk=0; kk<2; kk++) 
        ny += (nodal_phi[i+ii][j+1][k+kk] - nodal_phi[i+ii][j][k+kk]); 
    ny *= (0.25 / phi.dyc(j));

    real nz = 0.0;
    for(int ii=0; ii<2; ii++) 
      for(int jj=0; jj<2; jj++) 
        nz += (nodal_phi[i+ii][j+jj][k+1] - nodal_phi[i+ii][j+jj][k]); 
    nz *= (0.25 / phi.dzc(k));

    /* cell coordinates (node coordinates) */
    real x[2] = {phi.xn(i), phi.xn(i+1)};
    real y[2] = {phi.yn(j), phi.yn(j+1)};
    real z[2] = {phi.zn(k), phi.zn(k+1)};
 
    assert( phi.xc(i) > x[0] && phi.xc(i) < x[1] );

    CutCell * ccell = NULL;
    int n=0;

    /*------------------------+ 
    |  cuts in "i" direction  |
    +------------------------*/
    for(int jj=0; jj<2; jj++)
      for(int kk=0; kk<2; kk++) {

        real phi0 = nodal_phi[i]  [j+jj][k+kk];
        real phi1 = nodal_phi[i+1][j+jj][k+kk];
      
        real xi;
        if( phi0 < 0.5 && phi1 > 0.5 || phi0 > 0.5 && phi1 < 0.5) {
          xi = x[0] + (x[1] - x[0]) * (0.5 - phi0) / (phi1 - phi0);
          assert( xi > x[0] && xi < x[1] );
          xp[n] = xi; yp[n] = phi.yn(j+jj); zp[n] = phi.zn(k+kk);
          n++;
        }
      }

    /*------------------------+ 
    |  cuts in "j" direction  |
    +------------------------*/
    for(int ii=0; ii<2; ii++)
      for(int kk=0; kk<2; kk++) {

        real phi0 = nodal_phi[i+ii][j]  [k+kk];
        real phi1 = nodal_phi[i+ii][j+1][k+kk];
      
        real yi;
        if( phi0 < 0.5 && phi1 > 0.5 || phi0 > 0.5 && phi1 < 0.5) {
          yi = y[0] + (y[1] - y[0]) * (0.5 - phi0) / (phi1 - phi0);
          assert( yi > y[0] && yi < y[1] );
          xp[n] = phi.xn(i+ii); yp[n] = yi; zp[n] = phi.zn(k+kk);
          n++;
        }
      }

    /*------------------------+ 
    |  cuts in "k" direction  |
    +------------------------*/
    for(int ii=0; ii<2; ii++)
      for(int jj=0; jj<2; jj++) {

        real phi0 = nodal_phi[i+ii][j+jj][k];
        real phi1 = nodal_phi[i+ii][j+jj][k+1];
      
        real zi;
        if( phi0 < 0.5 && phi1 > 0.5 || phi0 > 0.5 && phi1 < 0.5) {
          zi = z[0] + (z[1] - z[0]) * (0.5 - phi0) / (phi1 - phi0);
          assert( zi > z[0] && zi < z[1] );
          xp[n] = phi.xn(i+ii); yp[n] = phi.yn(j+jj); zp[n] = zi;
          n++;
        }
      }

    /*-----------------------------------------------+ 
    |  if number of cutting nodes is between 3 and,  | 
    |  6 a regular cutting polygon can be created    |
    +-----------------------------------------------*/
    if(n >= 3 && n <= 6) {
      ccell = new CutCell();
      ccell->ijk(i,j,k);

      arrange_poly(n, xp,yp,zp, nx, ny, nz); 
      Polygon * cut_1 = new Polygon(n, xp, yp, zp); 

      polys.push_back(*cut_1);
    }

  }
}
 
/******************************************************************************/
void Body::cut(const Domain & dom) {

#ifdef DEBUG
  std::cout<<"body_cut::start. irank= "<<boil::cart.iam()<<"\n";
#endif
 
  /*-------------------+
  |  initialize array  |
  +-------------------*/
  cut_init(dom);
#ifdef DEBUG
  std::cout<<"body_cut::pass cut_init. irank= "<<boil::cart.iam()<<"\n";
#endif

  /*-----------------------------------------------------+
  |  cut staggered cells & calculate fv with flood fill  |
  +-----------------------------------------------------*/
  cut_stgd(dom);
#ifdef DEBUG
  std::cout<<"body_cut::pass cut_stgd. irank= "<<boil::cart.iam()<<"\n";
#endif

  /*--------------------------------------------------+
  |  cut scalar cells, cal. fv with flood fill, and   |
  |  replace original stl body with scalar cut faces  |
  +--------------------------------------------------*/
  cut_center(dom);
#ifdef DEBUG
  std::cout<<"body_cut::pass cut_center. irank= "<<boil::cart.iam()<<"\n";
#endif
  set_vecoff();

  /*------------------------------+
  |  calculate distance function  |
  +------------------------------*/
  distfunc(dom);
  nwall(dom);

  /*--------------------------------+
  |  set acceptor for cell merging  |
  +--------------------------------*/
#if 0
  acceptor();
#endif

  /*------------------------------------+
  |  store boundary cells data for CIP  |
  +------------------------------------*/
  //bd_center(dom, sca);

  /*----------------------------+
  |  essentially for debugging  |
  +----------------------------*/
#if 0
  if(boil::plot) boil::plot->plot((*sca), "sca");
  if(boil::plot) boil::plot->plot((*vec), "vec");
  if(boil::plot) boil::plot->plot(*this, "body");
  if(boil::plot) boil::plot->plot((*bdist), "bdist");
#endif

  /* release memory */
  //delete nfu,nfv,nfw;
#ifdef DEBUG
  std::cout<<"body_cut::end. irank= "<<boil::cart.iam()<<"\n";
#endif

}

/*-----------------------------------------------------------------------------+
 '$Id: body_cut.cpp,v 1.30 2014/08/06 08:25:22 sato Exp $'/
+-----------------------------------------------------------------------------*/
