#include "pathline.h"

/***************************************************************************//**
*  advance particles
*  Procedure: 1. calculate uvw
*             2. calculate additional variables to be traced
*******************************************************************************/
void Pathline::init() {

  real utmp[np()],vtmp[np()],wtmp[np()];
  *utmp = -boil::unreal;
  *vtmp = -boil::unreal;
  *wtmp = -boil::unreal;

  /* 1. calculate uvw */
  for (int ip = 0; ip < np(); ip++){
    // current position
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    // velocity at current position
    // uvw->interpolate = -1e+24, if particle #ip is not in decomposed domain
    utmp[ip] = uvw->interpolate(Comp::u(),xold,yold,zold);
    vtmp[ip] = uvw->interpolate(Comp::v(),xold,yold,zold);
    wtmp[ip] = uvw->interpolate(Comp::w(),xold,yold,zold);
    //if(ip==0){
    //  std::cout<<"advance: "<<boil::cart.iam()<<" "<<utmp[ip]<<"\n";
    //}
  }
  boil::cart.max_real_n(utmp,np());
  boil::cart.max_real_n(vtmp,np());
  boil::cart.max_real_n(wtmp,np());

  for (int ip = 0; ip < np(); ip++){
    // store velocity
    particles[ip].u(utmp[ip]);
    particles[ip].v(vtmp[ip]);
    particles[ip].w(wtmp[ip]);
    //if(ip==0){
    //  std::cout<<"advance: "<<boil::cart.iam()<<" "<<utmp[ip]<<"\n";
    //}
  }

  /* 2. calculate additional variables to be traced */
  for (int ival = 0; ival < nval(); ival++) {
    *utmp = -boil::unreal;
    for (int ip = 0; ip < np(); ip++){
      real xnew = particles[ip].x();
      real ynew = particles[ip].y();
      real znew = particles[ip].z();
      if (ival == 0) {
        utmp[ip] = s1->interpolate(xnew,ynew,znew);
      } else if (ival == 1) {
        utmp[ip] = s2->interpolate(xnew,ynew,znew);
      } else if (ival == 2) {
        utmp[ip] = s3->interpolate(xnew,ynew,znew);
      }
    }

    boil::cart.max_real_n(utmp,np());

    for (int ip = 0; ip < np(); ip++){
      particles[ip].sval(ival,utmp[ip]);
      //if(ip==0)boil::oout<<"sval="<<particles[ip].sval(ival)<<"\n";
    }
  }

  return;
}


