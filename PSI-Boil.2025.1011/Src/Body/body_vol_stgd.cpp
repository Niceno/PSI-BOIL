#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::vol_stgd(const Domain & dom, const Body & nfx, const Body & nfy,
                       const Body & nfz) {

  boil::timer.start("flood fill staggered");

  for_m(m)
    for_avmijk((*vec),m,i,j,k)
      (*vec)[m][i][j][k]=-1.0;

  /* put the volume fractions of cut cells 
     (a temporary measure, it can't stay like this) */
  for_m(m){
    for(int cc=0; cc<cells[~m].size(); cc++) {
      int i, j, k;
      cells[~m][cc].ijk(&i, &j, &k);
      (*vec)[m][i][j][k] = cells[~m][cc].fV();
    }
  }

  /* preparation for body-fill */
  stack_pointer = 0;
  stack_size    = 4 * dom.ni() * dom.nj() * dom.nk();

  if (stack_size>4*500*500*500) {
    boil::aout<<"body_vol_stgd:lack of memory. Use more cores!\n";
    exit(0);
  }

  alloc2d( & stack, stack_size, 3 );

  for_m(m){
    int iundef = 1;      // number of undefined cell
    int nproc = boil::cart.nproc();
    int irank = boil::cart.iam();
    int ipcal[nproc];
    int npcal;
    for(int i=0; i<nproc; i++)ipcal[i]=0;
    if (nccells_in[~m]>0) {
      do {
        /* find seed */
        int iseed, jseed, kseed, icseed;
        iseed = jseed = kseed = icseed = -1;
        int i, j, k;
        for(int cc=0; cc<nccells_in[~m]; cc++) {
          cells[~m][cc].ijk(&i, &j, &k);
          /* check neighbour */
          if( i-1>=sifl[~m]){
            if( (*vec)[m][i-1][j][k]<-0.5 ){
              iseed=i-1; jseed=j; kseed=k; icseed=cc;
              goto findseed;
            }
          }
          if( i+1<=eifl[~m]){
            if( (*vec)[m][i+1][j][k]<-0.5 ){
              iseed=i+1; jseed=j; kseed=k; icseed=cc;
              goto findseed;
            }
          }
          if( j-1>=sjfl[~m]){
            if( (*vec)[m][i][j-1][k]<-0.5 ){
              iseed=i; jseed=j-1; kseed=k; icseed=cc;
              goto findseed;
            }
          }
          if( j+1<=ejfl[~m]){
            if( (*vec)[m][i][j+1][k]<-0.5 ){
              iseed=i; jseed=j+1; kseed=k; icseed=cc;
              goto findseed;
            }
          }
          if( k-1>=skfl[~m]){
            if( (*vec)[m][i][j][k-1]<-0.5 ){
              iseed=i; jseed=j; kseed=k-1; icseed=cc;
              goto findseed;
            }
          }
          if( k+1<=ekfl[~m]){
            if( (*vec)[m][i][j][k+1]<-0.5 ){
              iseed=i; jseed=j; kseed=k+1; icseed=cc;
              goto findseed;
            }
          }
        }
        std::cout<<"body_fill_stgd:Error! Could not find! "<<irank<<"\n";
        findseed:

        /* flood fill */
        assert( iseed!=-1 ); assert( jseed!=-1 ); assert( kseed!=-1 );

        real color=-1;
        // estimate solid or fluid
        real nxf, nyf, nzf;
        if(m==Comp::u()){
          nxf=nfx[icseed].n(0);
          nyf=nfx[icseed].n(1);
          nzf=nfx[icseed].n(2);
        } else if(m==Comp::v()){
          nxf=nfy[icseed].n(0);
          nyf=nfy[icseed].n(1);
          nzf=nfy[icseed].n(2);
        } else {
          nxf=nfz[icseed].n(0);
          nyf=nfz[icseed].n(1);
          nzf=nfz[icseed].n(2);
        }
        real nxc=(*vec).xc(m,iseed)-(*vec).xc(m,i);
        real nyc=(*vec).yc(m,jseed)-(*vec).yc(m,j);
        real nzc=(*vec).zc(m,kseed)-(*vec).zc(m,k);
        if( nxf*nxc + nyf*nyc + nzf*nzc >= 0){
          color=1.0;
        } else {
          color=0.0;
        }

        body_fillv(color, iseed, jseed, kseed, m);

        /* count iundef */
        iundef=0;
        for(int i=sifl[~m]; i<=eifl[~m]; i++){
          for(int j=sjfl[~m]; j<=ejfl[~m]; j++){
            for(int k=skfl[~m]; k<=ekfl[~m]; k++){
              if( (*vec)(m)[i][j][k]<-0.5 ) iundef++;
            }
          }
        }
#ifdef DEBUG
        std::cout<<"vol_stgd::iundef= "<<iundef<<" "<<irank<<" "<<m<<"\n";
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
    if( nproc==1 && npcal==0){
      boil::oout<<"body_vol_stgd: Error! No interface cell\n";
      boil::oout<<"nccells_in[~m]= "<<nccells_in[~m]<<" "<<m<<"\n";
      exit(0);
    }
#ifdef DEBUG
    std::cout<<"body_vol_stgd:npcal= "<<npcal<<" "<<nproc<<" "<<m<<"\n";
    std::cout<<"nccells_in[m]= "<<nccells_in[~m]<<"\n";
#endif

    /* some procs are not finished */
    int iloop = 0;
    while (npcal < nproc){
      iloop++;
#ifdef DEBUG
      boil::aout<<"body_vol_stgd:iloop= "<<iloop<<" "<<irank<<"\n";
#endif
      /* exchange data */
      (*vec).exchange(ipcal,m);

      if(ipcal[irank]==0){
        /* flood fill2 */
        real color=-1.0;
        int sia=(*vec).si(m);
        int sja=(*vec).sj(m);
        int ska=(*vec).sk(m);
        int eia=(*vec).ei(m);
        int eja=(*vec).ej(m);
        int eka=(*vec).ek(m);
        if((*vec)[m][sia-1][sja][ska]>-0.5) color = (*vec)[m][sia-1][sja][ska];
        if((*vec)[m][sia][sja-1][ska]>-0.5) color = (*vec)[m][sia][sja-1][ska];
        if((*vec)[m][sia][sja][ska-1]>-0.5) color = (*vec)[m][sia][sja][ska-1];
        if((*vec)[m][eia+1][eja][eka]>-0.5) color = (*vec)[m][eia+1][eja][eka];
        if((*vec)[m][eia][eja+1][eka]>-0.5) color = (*vec)[m][eia][eja+1][eka];
        if((*vec)[m][eia][eja][eka+1]>-0.5) color = (*vec)[m][eia][eja][eka+1];
        if( color > -0.5){
          for_vmijk((*vec),m,i,j,k){
            (*vec)[m][i][j][k]=color;
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

    (*vec)(m).bnd_update();
    (*vec).exchange(m);

  }

  dealloc2d( & stack );
  boil::timer.stop("flood fill staggered");

#ifdef DEBUG
  std::cout<<"body_vol_stgd:end\n";
#endif

}
