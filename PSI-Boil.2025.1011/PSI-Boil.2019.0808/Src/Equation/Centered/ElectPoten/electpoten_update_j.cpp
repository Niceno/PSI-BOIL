#include "electpoten.h"
/***************************************************************************//**
*  Calculate electric current J = -sigma_e.grad(pot) + sigma_e.u x B
*                               = sigma_e.(-grad(pot) + u x B)
*******************************************************************************/
void ElectPoten::update_j() {

  phi.bnd_update();

  //real ufc,vfc,wfc,bx,by,bz;
  //real vect_A[3],vect_B[3],vect_C[3];

  for_m(m){
    for_vmijk(J,m,i,j,k){
      real vect_A[3],vect_B[3],vect_C[3];
      vect_A[0]=uf[m][i][j][k];
      vect_A[1]=vf[m][i][j][k];
      vect_A[2]=wf[m][i][j][k];
      valueFace(&B,m,i,j,k,&vect_B[0],&vect_B[1],&vect_B[2]);
      boil::crossProduct(vect_C, vect_A, vect_B);
      J[m][i][j][k] = fluid()->sigma_e(m,i,j,k)  
                    * (-phi.grad_face(m,i,j,k) + vect_C[~m]);
                  //* (-phi.grad_face(m,i,j,k) + (vfc*bz - wfc*by));
      //if(m==Comp::w()&&i==5&&k==5&&(j==3||j==4||j==5||j==2)){
      //  std::cout<<"update_j: "<<j<<" "<<J[m][i][j][k]<<" "
                 //<<phi.grad_face(m,i,j,k)<<" "<<vect_C[~m]<<"\n";
      //           <<vect_A[0]<<" "<<vect_A[1]<<" "<<vect_A[2]<<" "
      //           <<vect_B[0]<<" "<<vect_B[1]<<" "<<vect_B[2]<<"\n";
      //}
    }
  }

#if 0
  Comp m;
  m = Comp::u();
  for_vmijk(J,m,i,j,k){
    ufc=uf[m][i][j][k];
    vfc=vf[m][i][j][k];
    wfc=wf[m][i][j][k];
    valueFace(&B,m,i,j,k,&bx,&by,&bz);
    J[m][i][j][k] = fluid()->sigma_e(m,i,j,k) 
                  * (-phi.grad_face(m,i,j,k) + (vfc*bz - wfc*by));
  }

  m = Comp::v();
  for_vmijk(J,m,i,j,k){
    ufc=uf[m][i][j][k];
    vfc=vf[m][i][j][k];
    wfc=wf[m][i][j][k];
    valueFace(&B,m,i,j,k,&bx,&by,&bz);
    J[m][i][j][k] = fluid()->sigma_e(m,i,j,k)  
                  * (-phi.grad_face(m,i,j,k) + (wfc*bx - ufc*bz));
  }

  m = Comp::w();
  for_vmijk(J,m,i,j,k){
    ufc=uf[m][i][j][k];
    vfc=vf[m][i][j][k];
    wfc=wf[m][i][j][k];
    valueFace(&B,m,i,j,k,&bx,&by,&bz);
    J[m][i][j][k] = fluid()->sigma_e(m,i,j,k)  
                  * (-phi.grad_face(m,i,j,k) + (ufc*by - vfc*bx));
  }
#endif

  J.bnd_update_nooutlet();
  J.exchange_all();
  //std::cout<<"update_j:end\n";

  return;
}
