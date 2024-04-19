#include "pathline.h"

/***************************************************************************//**
*  advance particles
*******************************************************************************/
void Pathline::advance() {

  real dt = time->dt();
  for (int ip = 0; ip < np(); ip++){
    //std::cout<<"Pathline:advance:ip= "<<ip<<" x= "<<particles[ip].x()<<" y= "
    //       <<particles[ip].y()<<" z= "<<particles[ip].z()<<"\n";

    // current position
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();

    if ( uvw->domain()->global_min_x() > xold
      || uvw->domain()->global_max_x() < xold
      || uvw->domain()->global_min_y() > yold
      || uvw->domain()->global_max_y() < yold
      || uvw->domain()->global_min_z() > zold
      || uvw->domain()->global_max_z() < zold ) {
      boil::oout<<"pathline_advance:WARNING!!! Particle(no."<<ip
                <<") is outside of the domain. x,y,z= "
                <<xold<<" "<<yold<<" "<<zold<<"\n";
      particles[ip].u(0.0);
      particles[ip].v(0.0);
      particles[ip].w(0.0);

    } else {

      // velocity at current position
      real uu = uvw->Interpolate(Comp::u(),xold,yold,zold);
      real vv = uvw->Interpolate(Comp::v(),xold,yold,zold);
      real ww = uvw->Interpolate(Comp::w(),xold,yold,zold);

      // new position
      real xnew = xold + uu*dt;
      real ynew = yold + vv*dt;
      real znew = zold + ww*dt;
      particles[ip].x(xnew);
      particles[ip].y(ynew);
      particles[ip].z(znew);

      // store velocity
      particles[ip].u(uu);
      particles[ip].v(vv);
      particles[ip].w(ww);

      //boil::oout<<"pathline_advance:ip="<<ip<<" "<<time->current_time()
      //          <<" "<<xnew<<" "<<ynew<<" "<<znew<<"\n";
    }
  }
  return;
}


