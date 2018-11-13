#include "vof.h"

real VOF::calc_alpha(real & v, real & vma, real & vmb, real & vmc){

  real w, v1, v3, vm1, vm2, vm3, vm12; 
  real a0, a1, a2, q0, sp, th;
  real alpha;


  w     = boil::minr(v, 1.0-v);  
  vm1   = boil::minr(vma, vmb, vmc);
  vm3   = boil::maxr(vma, vmb, vmc);
  vm2   = fabs(1.0 - vm3 - vm1);
  vm12  = vm1 + vm2;
//  std::cout<<"w "<<w<<"\n";
  #if 0
  std::cout<<"w "<<w<<"\n";
  std::cout<<"vm1 "<<vm1<<"\n";
  std::cout<<"vm2 "<<vm2<<"\n";
  std::cout<<"vm3 "<<vm3<<"\n";
  std::cout<<"vm12 "<<vm12<<"\n";
  #endif
  v1    = pow(vm1, 2.0) / (6.0 * vm2 * vm3 + boil::pico);
  #if 0
  std::cout<<"v1 "<<v1<<"\n";
  std::cout<<"v2 "<<v1 + (vm2-vm1) / (2.0 * vm3)<<"\n";
  std::cout<<"v3 "<< 0.5 * vm12 / vm3 <<"\n";
  #endif

  if (w < v1){
    alpha = pow(6.0 * vm1 * vm2 * vm3 * w, 1.0/3.0);
  } else if(w < v1 + (vm2-vm1) / (2.0 * vm3)){
    alpha = 0.5 * (vm1 + sqrt(pow(vm1, 2) + 8.0 * vm2 * vm3 * (w - v1)));
  } else {
    alpha = 0.0;
    if (vm3 < vm12){
      v3 = (pow(vm3, 2.0) * (3.0 * vm12 - vm3) + pow(vm1, 2.0) * (vm1 - 3.0 *
           vm3) + pow(vm2, 2.0) * (vm2 - 3.0 * vm3)) / (6.0 * vm1 * vm2 * vm3);
    } else{
      v3 = 0.5 * vm12 / vm3;
      if(v3 <w || approx(v3,w)){
        alpha = vm3 * w + 0.5 * vm12;
      }
    }
//  std::cout<<"alpha "<< alpha  <<"\n";
    if (alpha == 0.0){
      if(w < v3){
        a2 = -3.0 * vm12;
        a1 = 3.0 * (pow(vm1, 2.0) + pow(vm2, 2.0));               
        a0 = -(pow(vm1, 3.0) + pow(vm2, 3.0) - 6.0 * vm1 * vm2 * vm3 * w);
      } else {
        a2 = -1.5;
        a1 = 1.5 * (pow(vm1, 2.0) + pow(vm2, 2.0) + pow(vm3, 2.0));
        a0 = -0.5 * (pow(vm1, 3.0) + pow(vm2, 3.0) + pow(vm3, 3.0) - 6.0 *
             vm1 * vm2 * vm3 * w);
      }
      q0     = (1.0/6) * (a1 * a2 - 3.0 * a0) - pow(a2, 3.0) * (1.0/27.0);
      sp    = sqrt(-1.0/3.0 * a1 + 1.0/9.0 * pow(a2, 2.0));
      th    = 1.0/3.0 * acos(q0 / (pow(sp, 3.0)));
      alpha = 2.0 * sp * cos(th + (4.0/3.0 * (boil::pi))) - (1.0/3.0) * a2;
    }
  }
 

  if(v > 0.5){
      alpha = 1.0 - alpha;
  }
  
  return alpha;
}

