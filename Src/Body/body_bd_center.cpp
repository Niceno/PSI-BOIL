#include "body.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"
#define DEBUG

/******************************************************************************/
void Body::bd_center(const Domain & dom, Scalar * sca) {
/***************************************************************************//**
*  \brief  set boundary-cell data which will be used in CIPCSL2
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"bd_center: begin.\n";
#endif
  int nccells=0;
  int i,j,k;

  /* i-min */
  if((dom.dim(Comp::i()) == 1 && dom.period(0)!=1) ||
     (dom.dim(Comp::i()) >= 1 && dom.neighbour(Dir::imin()) == par_proc_null)
     && dom.period(0)!=1){
    i=(*sca).si();
    for_vj( (*sca), j) {
      for_vk( (*sca), k) {
        if( (*sca)[i][j][k]>0.0 && (*sca)[i][j][k]<1.0){
          int ic=index[3][i][j][k];
          CutCell *cc = new CutCell();
          // set fP
          for(int jj=0; jj<=1; jj++){
            for(int kk=0; kk<=1; kk++){
              cc->fP(cells[3][ic].fP(0,jj,kk),1,jj,kk);
            }
          }
          // set fE
          cc->fE(cells[3][ic].fE(1,0,0),1,1,0);
          cc->fE(cells[3][ic].fE(1,0,1),1,1,1);
          cc->fE(cells[3][ic].fE(2,0,0),2,1,0);
          cc->fE(cells[3][ic].fE(2,0,1),2,1,1);
          // set fS
          cc->fS(cells[3][ic].fS(0,0),0,1);
          // store
          cc->ijk(i-1,j,k); cells[3].push_back( *cc );
          nccells++;
        }
      }
    }
  }

#ifdef DEBUG
  std::cout<<"bd_center: i-max\n";
#endif
  /* i-max */
  if((dom.dim(Comp::i()) == 1 && dom.period(0)!=1) ||
    ((dom.dim(Comp::i()) >= 1 && dom.neighbour(Dir::imax()) == par_proc_null)
     && dom.period(0)!=1)){
    i=(*sca).ei();
    for_vj( (*sca), j) {
      for_vk( (*sca), k) {
        if( (*sca)[i][j][k]>0.0 && (*sca)[i][j][k]<1.0){
          int ic=index[3][i][j][k];
          CutCell *cc = new CutCell();
          // set fP
          for(int jj=0; jj<=1; jj++){
            for(int kk=0; kk<=1; kk++){
              cc->fP(cells[3][ic].fP(1,jj,kk),0,jj,kk);
            }
          }
          // set fE
          cc->fE(cells[3][ic].fE(1,1,0),1,0,0);
          cc->fE(cells[3][ic].fE(1,1,1),1,0,1);
          cc->fE(cells[3][ic].fE(2,1,0),2,0,0);
          cc->fE(cells[3][ic].fE(2,1,1),2,0,1);
          // set fS
          cc->fS(cells[3][ic].fS(0,1),0,0);
          // store
          cc->ijk(i+1,j,k); cells[3].push_back( *cc );
          nccells++;
        }
      }
    }
  }

  for_avijk((*sca),i,j,k) index[3][i][j][k] = -1;
  for(int cc=0; cc<cells[3].size(); cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    index[3][i][j][k] = cc;
  }

#ifdef DEBUG
  std::cout<<"bd_center: i-max\n";
#endif
  /* j-min */
  if((dom.dim(Comp::j())==1 && dom.period(1)!=1) ||
    ((dom.dim(Comp::j()) >= 1 && dom.neighbour(Dir::jmin()) == par_proc_null)
     && dom.period(1)!=1)){
    std::cout<<"bd:jmin \n";
    j=(*sca).sj();
    for_avi( (*sca), i) {
      for_vk( (*sca), k) {
        if( (*sca)[i][j][k]>0.0 && (*sca)[i][j][k]<1.0){
          int ic=index[3][i][j][k];
          CutCell *cc = new CutCell();
          // set fP
          for(int ii=0; ii<=1; ii++){
            for(int kk=0; kk<=1; kk++){
              cc->fP(cells[3][ic].fP(ii,0,kk),ii,1,kk);
            }
          }
          // set fE
          cc->fE(cells[3][ic].fE(0,0,0),0,1,0);
          cc->fE(cells[3][ic].fE(0,0,1),0,1,1);
          cc->fE(cells[3][ic].fE(2,0,0),2,0,1);
          cc->fE(cells[3][ic].fE(2,1,0),2,1,1);
          // set fS
          cc->fS(cells[3][ic].fS(1,0),1,1);
          // store
          cc->ijk(i,j-1,k); cells[3].push_back( *cc );
          nccells++;
        }
      }
    }
  }

#ifdef DEBUG
  std::cout<<"bd_center: j-max\n";
#endif
  /* j-max */
  if((dom.dim(Comp::j()) == 1 && dom.period(1)!=1) ||
    ((dom.dim(Comp::j()) >= 1 && dom.neighbour(Dir::jmax()) == par_proc_null)
     && dom.period(1)!=1)){
    j=(*sca).ej();
    for_avi( (*sca), i) {
      for_vk( (*sca), k) {
        if( (*sca)[i][j][k]>0.0 && (*sca)[i][j][k]<1.0){
          int ic=index[3][i][j][k];
          CutCell *cc = new CutCell();
          // set fP
          for(int ii=0; ii<=1; ii++){
            for(int kk=0; kk<=1; kk++){
              cc->fP(cells[3][ic].fP(ii,1,kk),ii,0,kk);
            }
          }
          // set fE
          cc->fE(cells[3][ic].fE(0,1,0),0,0,0);
          cc->fE(cells[3][ic].fE(0,1,1),0,0,1);
          cc->fE(cells[3][ic].fE(2,0,1),2,0,0);
          cc->fE(cells[3][ic].fE(2,1,1),2,1,0);
          // set fS
          cc->fS(cells[3][ic].fS(1,1),1,0);
          // store
          cc->ijk(i,j+1,k); cells[3].push_back( *cc );
          nccells++;
        }
      }
    }
  }

  for_avijk((*sca),i,j,k) index[3][i][j][k] = -1;
  for(int cc=0; cc<cells[3].size(); cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    index[3][i][j][k] = cc;
  }

#ifdef DEBUG
  std::cout<<"bd_center: k-min\n";
#endif
  /* k-min */
  if((dom.dim(Comp::k()) == 1 && dom.period(2)!=1) ||
    ((dom.dim(Comp::k()) >= 1 && dom.neighbour(Dir::kmin()) == par_proc_null)
     && dom.period(2)!=1)){
    k=(*sca).sk();
    std::cout<<"sk= "<<k<<"\n";
    for_avi( (*sca), i) {
      for_avj( (*sca), j) {
        if( (*sca)[i][j][k]>0.0 && (*sca)[i][j][k]<1.0){
          int ic=index[3][i][j][k];
          CutCell *cc = new CutCell();
          // set fP
          for(int ii=0; ii<=1; ii++){
            for(int jj=0; jj<=1; jj++){
              cc->fP(cells[3][ic].fP(ii,jj,0),ii,jj,1);
            }
          }
          // set fE
          cc->fE(cells[3][ic].fE(0,0,0),0,0,1);
          cc->fE(cells[3][ic].fE(0,1,0),0,1,1);
          cc->fE(cells[3][ic].fE(1,0,0),1,0,1);
          cc->fE(cells[3][ic].fE(1,1,0),1,1,1);
          // set fS
          cc->fS(cells[3][ic].fS(2,0),2,1);
          // store
          cc->ijk(i,j,k-1); cells[3].push_back( *cc );
          nccells++;
        }
      }
    }
  }

#ifdef DEBUG
  std::cout<<"bd_center: k-max\n";
#endif
  /* k-max */
  if((dom.dim(Comp::k()) == 1 && dom.period(2)!=1) ||
    ((dom.dim(Comp::k()) >= 1 && dom.neighbour(Dir::kmax()) == par_proc_null)
     && dom.period(2)!=1)){
    k=(*sca).ek();
    for_avi( (*sca), i) {
      for_avj( (*sca), j) {
        if( (*sca)[i][j][k]>0.0 && (*sca)[i][j][k]<1.0){
          int ic=index[3][i][j][k];
          CutCell *cc = new CutCell();
          // set fP
          for(int ii=0; ii<=1; ii++){
            for(int jj=0; jj<=1; jj++){
              cc->fP(cells[3][ic].fP(ii,jj,1),ii,jj,0);
            }
          }
          // set fE
          cc->fE(cells[3][ic].fE(0,0,1),0,0,0);
          cc->fE(cells[3][ic].fE(0,1,1),0,1,0);
          cc->fE(cells[3][ic].fE(1,0,1),1,0,0);
          cc->fE(cells[3][ic].fE(1,1,1),1,1,0);
          // set fS
          cc->fS(cells[3][ic].fS(2,1),2,0);
          // store
          cc->ijk(i,j,k+1); cells[3].push_back( *cc );
          nccells++;
        }
      }
    }
  }

  for_avijk((*sca),i,j,k) index[3][i][j][k] = -1;
  for(int cc=0; cc<cells[3].size(); cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    index[3][i][j][k] = cc;
  }

  nccells_bd[3] += nccells;

#ifdef DEBUG
  std::cout<<"body_cut,nccells_in,_bd,size= "<<nccells_in[3]<<" "
           <<nccells_bd[3]<<" "<<cells[3].size()<<"\n";
#endif

}

/*-----------------------------------------------------------------------------+
 '$Id: body_bd_center.cpp,v 1.4 2014/02/03 14:12:33 sato Exp $'/
+-----------------------------------------------------------------------------*/
