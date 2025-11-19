#include "vof.h"

real VOF::calc_v(real alpha, real vma, real vmb, real vmc){

  real v;
  real a;  
  real vm1, vm2, vm3, vm12;

  a   = boil::minr(alpha, 1.0-alpha);
  v   = 0.0; 

  if(a > 0){
    vm1  = boil::minr(vma, vmb, vmc);
    vm3  = boil::maxr(vma, vmb, vmc);
    vm2  = fabs(1.0 - vm3 - vm1);
    vm12 = vm1 + vm2;
    if(a < vm1){
      v = pow(a, 3.0) / (6.0 * vm1 * vm2 * vm3);
    } else if(a < vm2){
      v = a * (a - vm1) / (2.0 * vm2 * vm3) + pow(vm1, 2.0) / (6.0 * vm2 *
          vm3 + boil::pico);
    } else if(a < boil::minr(vm12, vm3)){
      v = (pow(a, 2.0) * (3.0 * vm12 - a) + pow(vm1, 2.0) *  (vm1 - 3.0 * a) +
           pow(vm2, 2.0) * (vm2 -3.0 * a)) / (6.0 * vm1 * vm2 * vm3);
    } else if(vm3 < vm12){
      v = (pow(a, 2.0) * (3.0 - 2.0 * a) + pow(vm1, 2.0) * (vm1 - 3.0 * a) +
          pow(vm2, 2.0) * (vm2 - 3.0 * a) + pow(vm3, 2.0) * (vm3 - 3.0 * a)) /
          (6.0 * vm1 * vm2 * vm3);
    } else{
      v = (a - 0.5 * vm12) / vm3;
    }
  }
  if(alpha > 0.5){
    v = 1.0 - v;
  }
  
  return v;
}
