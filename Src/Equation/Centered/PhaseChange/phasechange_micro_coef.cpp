#include "phasechange.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
real PhaseChange::pcrate(real tpr) {
  real rate = 7.031e-5;   // A=3.6e-20, angle=80, y1=1.1e-3, y3=500
  real tpr2=tpr-100.0; 
  rate = rate / 9.0 * tpr2;
  //std::cout<<rate<<"\n";
  return 4.0*rate;
}

/******************************************************************************/
real PhaseChange::pcfrc(real tpr) {
  //real frc=0.0;
  real frc = - nucl->sigma() * (1.0 - cos(nucl->cangle()/180.0*pi));
  //std::cout<<"pcfrc: "<<frc<<"\n";
  return frc;
}

