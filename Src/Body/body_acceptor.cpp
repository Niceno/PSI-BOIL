#include "body.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"
#include "../Domain/domain.h"
#include <iomanip>
#define DEBUG

/******************************************************************************/
void Body::acceptor() {
/***************************************************************************//**
*  /brief Convert color function to distance function.
*         Reference:G.,Russo and P.,Smereka,"A remark on computing distance
*         functions",J.Comp.Phys.,Vol.163,2000,pp.51-67
*           input    :sca
*           output   :bdist
*           temporary:stmp,dflag
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"body_acceptor::start "<<boil::cart.iam()<<"\n";
#endif

  for(int cc=0; cc<nccells_in[3]; cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    if (cells[3][cc].fV()<0.5 && cells[3][cc].fV()>0.0) { // donor cell
      //std::cout<<i<<" "<<j<<" "<<k<<" "<<cells[3][cc].fV()<<"\n";
      real nx = ((*bdist)[i+1][j][k]-(*bdist)[i-1][j][k]) 
              / (bdist->dxw(i)+bdist->dxe(i));
      real ny = ((*bdist)[i][j+1][k]-(*bdist)[i][j-1][k]) 
              / (bdist->dys(j)+bdist->dyn(j));
      real nz = ((*bdist)[i][j][k+1]-(*bdist)[i][j][k-1]) 
              / (bdist->dzb(k)+bdist->dzt(k));
      real anx = fabs(nx);
      real any = fabs(ny);
      real anz = fabs(nz);
      if (anx>=any && anx>=anz) {
        if(nx>0){
          cells[3][cc].set_iacpt(i+1);
        } else {
          cells[3][cc].set_iacpt(i-1);
        }
        cells[3][cc].set_jacpt(j);
        cells[3][cc].set_kacpt(k);
      } else if (any>=anx && any>=anz) {
        if(ny>0){
          cells[3][cc].set_jacpt(j+1); 
        } else {
          cells[3][cc].set_jacpt(j-1);
        }
        cells[3][cc].set_iacpt(i);
        cells[3][cc].set_kacpt(k);
      } else if (anz>=anx && anz>=any) {
        if(nz>0){
          cells[3][cc].set_kacpt(k+1); 
        } else {
          cells[3][cc].set_kacpt(k-1);
        }
        cells[3][cc].set_iacpt(i);
        cells[3][cc].set_jacpt(j);
      } else {
        std::cout<<"body_acceptor:Error! \n";
        exit(0);
      }
      std::cout<<"donor= "<<i<<" "<<j<<" "<<k<<" acceptor= "
               <<cells[3][cc].iacpt()<<" "
               <<cells[3][cc].jacpt()<<" "
               <<cells[3][cc].kacpt()<<"\n";
    }
  }

  /* check whether acceptor is donor */
#if 0  // this algorithm doe not work for domain decomposition
  for(int cc=0; cc<nccells_in[3]; cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    if (cells[3][cc].fV()<0.5) { // donor cell
      int n_error=0;
      do {
        int n_error=0;
        int iac = cells[3][cc].iacpt();
        int jac = cells[3][cc].jacpt();
        int kac = cells[3][cc].kacpt();
        int cc_new=index[3][iac][jac][kac];

        //if(cells[3][cc_new].fV()<0.5){  //acceptor = donor
        if((*sca)[iac][jac][kac]<0.5){  //acceptor = donor , consider MPI
          n_error=1;
          int iac_new = cells[3][cc_new].iacpt();
          int jac_new = cells[3][cc_new].jacpt();
          int kac_new = cells[3][cc_new].kacpt();
          cells[3][cc].set_iacpt(iac_new);
          cells[3][cc].set_jacpt(jac_new);
          cells[3][cc].set_kacpt(kac_new);
        }
      } while(n_error>=1);
    }
  }
#endif
}

/*-----------------------------------------------------------------------------+
 '$Id: body_acceptor.cpp,v 1.1 2011/03/28 08:12:04 sato Exp $'/
+-----------------------------------------------------------------------------*/
