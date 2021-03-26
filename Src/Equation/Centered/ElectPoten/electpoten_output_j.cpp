#include "electpoten.h"
/***************************************************************************//**
*  Output electric current
*    j_potential: -sigma_e.grad(pot)
*    j_lorenz:  sigma_e.u x B
*******************************************************************************/
void ElectPoten::output_j() {

  phi.bnd_update();

  for_m(m){
    for_vmijk(J,m,i,j,k){
      J[m][i][j][k] = -fluid()->sigma_e(m,i,j,k)*phi.grad_face(m,i,j,k);
    }
  }
  J.bnd_update_nooutlet();
  J.exchange_all();
  boil::plot->plot(J,phi, "Jpotential-pot", time->current_step());

  for_m(m){
    for_vmijk(J,m,i,j,k){
      real vect_A[3],vect_B[3],vect_C[3];
      vect_A[0]=uf[m][i][j][k];
      vect_A[1]=vf[m][i][j][k];
      vect_A[2]=wf[m][i][j][k];
      valueFace(&B,m,i,j,k,&vect_B[0],&vect_B[1],&vect_B[2]);
      boil::crossProduct(vect_C, vect_A, vect_B);
      J[m][i][j][k] = fluid()->sigma_e(m,i,j,k)*vect_C[~m];
    }
  }
  J.bnd_update_nooutlet();
  J.exchange_all();
  boil::plot->plot(J,phi, "Jlorenz-pot", time->current_step());


  // recalculate j
  update_j();

  return;
}
