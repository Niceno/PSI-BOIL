#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::vol_center(const Domain & dom, const Body & nf) {

  boil::timer.start("flood fill center");

  for_avijk((*sca),i,j,k)
    (*sca)[i][j][k] = -1.0;
    
  /* put the volume fractions of cut cells 
     (a temporary measure, it can't stay like this) */
  for(int cc=0; cc<cells[3].size(); cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    (*sca)[i][j][k] = cells[3][cc].fV();
  }
  /* preparation for body-fill */
  stack_pointer = 0;
  stack_size    = 4 * dom.ni() * dom.nj() * dom.nk();
  boil::oout<<"body_vol_center:alloc2d\n";
  alloc2d( & stack, stack_size, 3 );

  int iundef = 1;      // number of undefined cell
  int nproc = boil::cart.nproc();
  int irank = boil::cart.iam();
  int ipcal[nproc];
  int npcal;
  for(int i=0; i<nproc; i++)ipcal[i]=0;

  if (cells[3].size()>0) {
    do {
      /* find seed */
      int iseed, jseed, kseed, icseed;
      iseed = jseed = kseed = icseed = -1;
      int i, j, k;
      for(int cc=0; cc<cells[3].size(); cc++) {
        cells[3][cc].ijk(&i, &j, &k);
        /* check neighbour */
        if( i-1>=sifl[3]){
          if( (*sca)[i-1][j][k]<-0.5 ){
            iseed=i-1; jseed=j; kseed=k; icseed=cc;
            goto findseed;
          }
        }
        if( i+1<=eifl[3]){
          if( (*sca)[i+1][j][k]<-0.5 ){
            iseed=i+1; jseed=j; kseed=k; icseed=cc;
            goto findseed;
          }
        }
        if( j-1>=sjfl[3]){
          if( (*sca)[i][j-1][k]<-0.5 ){
            iseed=i; jseed=j-1; kseed=k; icseed=cc;
            goto findseed;
          }
        }
        if( j+1<=ejfl[3]){
          if( (*sca)[i][j+1][k]<-0.5 ){
            iseed=i; jseed=j+1; kseed=k; icseed=cc;
            goto findseed;
          }
        }
        if( k-1>=skfl[3]){
          if( (*sca)[i][j][k-1]<-0.5 ){
            iseed=i; jseed=j; kseed=k-1; icseed=cc;
            goto findseed;
          }
        }
        if( k+1<=ekfl[3]){
          if( (*sca)[i][j][k+1]<-0.5 ){
            iseed=i; jseed=j; kseed=k+1; icseed=cc;
            goto findseed;
          }
        }
      }
      std::cout<<"body_fill_center:Error! Could not find! "<<irank<<"\n";
      findseed:
      /* flood fill */
      assert( iseed!=-1 ); assert( jseed!=-1 ); assert( kseed!=-1 );

      real color=-1;
      // estimate solid or fluid
      real nxf=nf[icseed].n(0);
      real nyf=nf[icseed].n(1);
      real nzf=nf[icseed].n(2);
      real nxc=(*sca).xc(iseed)-(*sca).xc(i);
      real nyc=(*sca).yc(jseed)-(*sca).yc(j);
      real nzc=(*sca).zc(kseed)-(*sca).zc(k);
      if( nxf*nxc + nyf*nyc + nzf*nzc >= 0){
        color=1.0;
      } else {
        color=0.0;
      }

      body_fills(color, iseed, jseed, kseed);

      /* count iundef */
      iundef=0;
      for(int i=sifl[3]; i<=eifl[3]; i++){
        for(int j=sjfl[3]; j<=ejfl[3]; j++){
          for(int k=skfl[3]; k<=ekfl[3]; k++){
            if( (*sca)[i][j][k]<-0.5 ) iundef++;
          }
        }
      }
#if DEBUG
      std::cout<<"center_fill::iundef= "<<iundef<<" irank="<<irank<<"\n";
#endif

    } while(iundef>0);

    /* cal npcal */
    ipcal[irank]=1;

  }

  boil::cart.sum_int_n(ipcal,nproc);
  npcal=0;
  for(int i=0; i<nproc; i++){
    npcal += ipcal[i];
  }

  /* some procs are not finished */
  int iloop = 0;
  while (npcal < nproc){
    iloop++;

#if DEBUG
    boil::aout<<"body_vol_center:iloop= "<<iloop<<" "<<irank<<"\n";
#endif

    /* exchange data */
    (*sca).exchange(ipcal);

    if(ipcal[irank]==0){
      /* flood fill2 */
      real color=-1.0;
      int sia=(*sca).si();
      int sja=(*sca).sj();
      int ska=(*sca).sk();
      int eia=(*sca).ei();
      int eja=(*sca).ej();
      int eka=(*sca).ek();
      if( (*sca)[sia-1][sja][ska]>-0.5 ) color = (*sca)[sia-1][sja][ska];
      if( (*sca)[sia][sja-1][ska]>-0.5 ) color = (*sca)[sia][sja-1][ska];
      if( (*sca)[sia][sja][ska-1]>-0.5 ) color = (*sca)[sia][sja][ska-1];
      if( (*sca)[eia+1][eja][eka]>-0.5 ) color = (*sca)[eia+1][eja][eka];
      if( (*sca)[eia][eja+1][eka]>-0.5 ) color = (*sca)[eia][eja+1][eka];
      if( (*sca)[eia][eja][eka+1]>-0.5 ) color = (*sca)[eia][eja][eka+1];
      if( color > -0.5){
        for_vijk((*sca),i,j,k){
          (*sca)[i][j][k]=color;
        }
        ipcal[irank]=1;
      }
    }

    /* cal npcal */
    boil::cart.max_int_n(ipcal,nproc);
    npcal=0;
    for(int i=0; i<nproc; i++){
      npcal += ipcal[i];
    }
  }

  sca->bnd_update();
  sca->exchange_all();

  dealloc2d( & stack );

  boil::timer.stop("flood fill center");
}
