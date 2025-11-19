#include "electpoten.h"
/***************************************************************************//**
*  Calculate velocity at cell faace
*    input: (*u)
*    output: vf, vf,wf 
*  objectives: to impose exact boundary condition to velocity field 
*  which is used for the outer product
*
*  VELOCITY-U
*                |             |                        |             |
*  j         ----+------*------+----        k       ----+------*------+----
*  ^             |             |            ^           |             |
*  |             *<-uf[mi]     *            |           *<-uf[mi]     *
*  +-->i         |             |            +-->i       |             |
*            ----+------*------+----                ----+------*------+----
*                |   uf[mj]    |                        |    uf[mk]   |
*
*  VELOCITY-V
*                |             |                        |             |
*  j         ----+------*------+----        k       ----+------*------+----
*  ^             |             |            ^           |             |
*  |             *<-vf[mi]     *            |           *<-vf[mi]     *
*  +-->i         |             |            +-->i       |             |
*            ----+------*------+----                ----+------*------+----
*                |   vf[mj]    |                        |    vf[mk]   |
*
*  VELOCITY-W
*                |             |                        |             |
*  k         ----+------*------+----        k       ----+------*------+----
*  ^             |             |            ^           |             |
*  |             *<-wf[mi]     *            |           *<-wf[mj]     *
*  +-->i         |             |            +-->j       |             |
*            ----+------*------+----                ----+------*------+----
*                |   wf[mk]    |                        |    wf[mk]   |
*
*******************************************************************************/
void ElectPoten::update_vel() {

  Comp mi = Comp::u();
  Comp mj = Comp::v();
  Comp mk = Comp::w();

  /* value at i-face */
  //std::cout<<"range:si.mi= "<<vf.si(mi)<<" "<<(*u).si(mi)<<"\n";
  //std::cout<<"range:sj.mi= "<<vf.sj(mi)<<" "<<(*u).sj(mi)<<"\n";
  //std::cout<<"range:sk.mi= "<<vf.sk(mi)<<" "<<(*u).sk(mi)<<"\n";
  for_vmijk(vf,mi,i,j,k){
    uf[mi][i][j][k] = (*u)[mi][i][j][k];
  }
  for_vmijk(vf,mi,i,j,k){
    vf[mi][i][j][k] = 0.25 * ((*u)[mj][i-1][j][k]+(*u)[mj][i-1][j+1][k]
                             +(*u)[mj][i  ][j][k]+(*u)[mj][i  ][j+1][k]);
  }
  for_vmijk(wf,mi,i,j,k){
    wf[mi][i][j][k] = 0.25 * ((*u)[mk][i-1][j][k]+(*u)[mk][i-1][j][k+1]
                             +(*u)[mk][i  ][j][k]+(*u)[mk][i  ][j][k+1]);
  }

  /* value at j-face */
  for_vmijk(vf,mj,i,j,k){
    uf[mj][i][j][k] = 0.25 * ((*u)[mi][i][j-1][k]+(*u)[mi][i+1][j-1][k]
                             +(*u)[mi][i][j  ][k]+(*u)[mi][i+1][j  ][k]);
  }
  for_vmijk(vf,mj,i,j,k){
    vf[mj][i][j][k] = (*u)[mj][i][j][k];
  }
  for_vmijk(wf,mj,i,j,k){
    wf[mj][i][j][k] = 0.25 * ((*u)[mk][i][j-1][k]+(*u)[mk][i][j-1][k+1]
                             +(*u)[mk][i][j  ][k]+(*u)[mk][i][j  ][k+1]);
  }

  /* value at k-face */
  //std::cout<<"range:si.mk= "<<vf.si(mk)<<" "<<(*u).si(mk)<<"\n";
  //std::cout<<"range:sj.mk= "<<vf.sj(mk)<<" "<<(*u).sj(mk)<<"\n";
  //std::cout<<"range:sk.mk= "<<vf.sk(mk)<<" "<<(*u).sk(mk)<<"\n";
  for_vmijk(vf,mk,i,j,k){
    uf[mk][i][j][k] = 0.25 * ((*u)[mi][i][j][k-1]+(*u)[mi][i+1][j][k-1]
                             +(*u)[mi][i][j][k  ]+(*u)[mi][i+1][j][k  ]);
  }
  for_vmijk(vf,mk,i,j,k){
    vf[mk][i][j][k] = 0.25 * ((*u)[mj][i][j][k-1]+(*u)[mj][i][j+1][k-1]
                             +(*u)[mj][i][j][k  ]+(*u)[mj][i][j+1][k  ]);
  }
  for_vmijk(wf,mk,i,j,k){
    wf[mk][i][j][k] = (*u)[mk][i][j][k];
  }

#if 0
  std::cout<<"update_vel:u  "<<(*u)[mi][5][5][2]<<" "<<(*u)[mi][5][5][3]<<" "<<(*u)[mi][5][5][4]<<"\n";
  std::cout<<"update_vel:uf "<<uf[mk][5][5][2]<<" "<<uf[mk][5][5][3]<<" "<<uf[mk][5][5][4]<<"\n";
  std::cout<<"update_vel:vf "<<vf[mk][5][5][3]<<" "<<vf[mk][5][5][4]<<"\n";
  std::cout<<"update_vel:wf "<<wf[mk][5][5][3]<<" "<<wf[mk][5][5][4]<<"\n";
#endif

  uf.bnd_update_nooutlet();
  uf.exchange_all();
  vf.bnd_update_nooutlet();
  vf.exchange_all();
  wf.bnd_update_nooutlet();
  wf.exchange_all();

#if 0
  std::cout<<"update_vel:uf "<<uf[mk][5][5][2]<<" "<<uf[mk][5][5][3]<<" "<<uf[mk][5][5][4]<<"\n";
  std::cout<<"update_vel:vf "<<vf[mk][5][5][3]<<" "<<vf[mk][5][5][4]<<"\n";
  std::cout<<"update_vel:wf "<<wf[mk][5][5][3]<<" "<<wf[mk][5][5][4]<<"\n";
  exit(0);
#endif

  return;
}
