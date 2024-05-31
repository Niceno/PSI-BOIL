#include "solidparticle.h"

/***************************************************************************//**
*  advance solid particles
*  Procedure: 1. calculate uvw
*             2. update position
*             3. calculate uvw at new position
*             4. calculate additional variables to be traced
*******************************************************************************/
void SolidParticle::advance() {

  boil::timer.start("solidparticle advance");

  real dt = time->dt();

  /* density of fluid */
  real dflu[np()], mu[np()];
#if 1
  //*dflu = -boil::unreal; // density   BUG
  //*mu = -boil::unreal;   // kinematic viscosity  BUG
  std::fill(dflu, dflu + np(), -boil::unreal);
  std::fill(mu,   mu + np(), -boil::unreal);
#else
  for (int ip = 0; ip < np(); ip++){
    dflu[ip]=-boil::unreal;
    mu  [ip]=-boil::unreal;
  }
#endif
  for (int ip = 0; ip < np(); ip++){
    // current position
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    int I = fluid()->domain()->I(xold);
    int J = fluid()->domain()->J(yold);
    int K = fluid()->domain()->K(zold);
    if (fluid()->domain()->contains_IJK(I,J,K)) {
      int local_i = fluid()->domain()->local_i(I);
      int local_j = fluid()->domain()->local_j(J);
      int local_k = fluid()->domain()->local_k(K);
      //if (boil::cart.iam()!=0)
      //std::cout<<"density: "<<boil::cart.iam()<<" "<<I<<" "<<J<<" "<<K<<" "<<local_i<<" "<<local_j<<" "<<local_k<<" "<<fluid()->rho(local_i,local_j,local_k)<<"\n";
      dflu[ip]=fluid()->rho(local_i,local_j,local_k);
      mu[ip]  =fluid()->mu (local_i,local_j,local_k);
      if(dflu[ip]>1e+10){
        std::cout<<"dflu= "<<dflu[ip]<<" ip "<<ip<<" "<<boil::cart.iam()<<"";
        exit(0);
      }
    }
  }
  //int ip = 605;
  //std::cout<<"ip= "<<ip<<" dflu "<<dflu[ip]<<" "<<mu[ip]<<" "<<boil::cart.iam()<<"\n";

  boil::cart.max_real_n(dflu,np());
  boil::cart.max_real_n(mu,np());

  //ip = 605;
  //boil::oout<<"after: ip= "<<ip<<" dflu "<<dflu[ip]<<" "<<mu[ip]<<"\n";

  /* velocity of fluid */
  real uflu[np()],vflu[np()],wflu[np()];
  //*uflu = -boil::unreal;
  //*vflu = -boil::unreal;
  //*wflu = -boil::unreal;
  std::fill(uflu, uflu + np(), -boil::unreal);
  std::fill(vflu, vflu + np(), -boil::unreal);
  std::fill(wflu, wflu + np(), -boil::unreal);

  for (int ip = 0; ip < np(); ip++){
    // current position
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    // velocity at current position
    // uvw->interpolate = -1e+24, if particle #ip is not in decomposed domain
    uflu[ip] = uvw->interpolate(Comp::u(),xold,yold,zold);
    vflu[ip] = uvw->interpolate(Comp::v(),xold,yold,zold);
    wflu[ip] = uvw->interpolate(Comp::w(),xold,yold,zold);
    //if(ip==0){
    //  std::cout<<"advance: "<<boil::cart.iam()<<" "<<utmp[ip]<<"\n";
    //}
  }
  boil::cart.max_real_n(uflu,np());
  boil::cart.max_real_n(vflu,np());
  boil::cart.max_real_n(wflu,np());

  /* equation of motion */
  for (int ip = 0; ip < np(); ip++){
    // mass
    real dia = particles[ip].diameter();
    real rad = 0.5 * dia;
    real vol = 4.0 / 3.0 * boil::pi * pow(rad,3.0);
    real mass = particles[ip].density() * vol;
    if(mass==0.0) {
      std::cout<<"msol=0.0 "<<ip<<" "<<boil::cart.iam()<<"\n";
      exit(0);
    }
    // Re
    real up = particles[ip].u();
    real vp = particles[ip].v();
    real wp = particles[ip].w();
    real velmag = sqrt(pow(uflu[ip]-up,2.0)
                     + pow(vflu[ip]-vp,2.0)
                     + pow(wflu[ip]-wp,2.0));
    real Re = dflu[ip] * dia * velmag / mu[ip]; 
    Re = boil::maxr(0.01, Re);
    // Cd
    real Cd = 0.44;
    if (Re < 1000) {
      Cd = 24.0 / Re * (1.0 + 0.15 * pow(Re,0.687));
    }
    // Drag force
    real A = rad * rad * boil::pi;
    real Fdx = A * 0.5 * Cd * dflu[ip] * velmag * ((uflu[ip]-up));
    real Fdy = A * 0.5 * Cd * dflu[ip] * velmag * ((vflu[ip]-vp));
    real Fdz = A * 0.5 * Cd * dflu[ip] * velmag * ((wflu[ip]-wp));
    // Draf force limit = m * (up - uf)/dt
    // Derivation:     F = m * a
    //            F * dt = m * a * dt
    //            F * dt = m * (u2 - u1), where u2 is the updated velocity
    //            The drag force makes u2 = u_fluid in the maximum condition
    real Fdx_limit = mass * (uflu[ip]-up)/dt;
    real Fdy_limit = mass * (vflu[ip]-vp)/dt;
    real Fdz_limit = mass * (wflu[ip]-wp)/dt;
    // Limit Fd
    Fdx = copysign(1.0,Fdx) * boil::minr(fabs(Fdx),fabs(Fdx_limit));
    Fdy = copysign(1.0,Fdy) * boil::minr(fabs(Fdy),fabs(Fdy_limit));
    Fdz = copysign(1.0,Fdz) * boil::minr(fabs(Fdz),fabs(Fdz_limit));
    // buoyancy force
    real den_diff = particles[ip].density() - dflu[ip];
    real Fbx = vol * gravity_x * den_diff;
    real Fby = vol * gravity_y * den_diff;
    real Fbz = vol * gravity_z * den_diff;
    // acceleration
    real accx = (Fdx + Fbx)/mass;
    real accy = (Fdy + Fby)/mass;
    real accz = (Fdz + Fbz)/mass;
    // velocity
    real up_new = up + dt * accx;
    real vp_new = vp + dt * accy;
    real wp_new = wp + dt * accz;
    // store velocity
    particles[ip].u(up_new);
    particles[ip].v(vp_new);
    particles[ip].w(wp_new);
    // current position
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    // new position
    real xnew = xold + dt * 0.5*(up+up_new);
    real ynew = yold + dt * 0.5*(vp+vp_new);
    real znew = zold + dt * 0.5*(wp+wp_new);
#if 0
    if(ip==547)
    boil::oout<<" ip "<<ip<<" mass "<<mass<<" Re "<<Re<<" Cd "<<Cd
              <<" up "<<up<<" vp "<<vp<<" wp "<<wp
              <<" uf "<<uflu[ip]<<" vf "<<vflu[ip]<<" wf "<<wflu[ip]
	      <<" df "<<dflu[ip]<<" mu "<<mu[ip]<<" velmag "<<velmag
	      <<" Fdx "<<Fdx<<" Fdy "<<Fdy<<" Fdz "<<Fdz
	      <<" Fbx "<<Fbx<<" Fby "<<Fby<<" Fbz "<<Fbz
	      <<" accx "<<accx<<" accy "<<accy<<" accz "<<accz
	      <<" dt "<<dt
	      <<" up_new "<<up_new<<" vp_new "<<vp_new<<" wp_new "<<wp_new
	      <<"\n";
#endif
    // check position inside the domain
    if ( uvw->domain()->global_min_x() > xnew
      || uvw->domain()->global_max_x() < xnew
      || uvw->domain()->global_min_y() > ynew
      || uvw->domain()->global_max_y() < ynew
      || uvw->domain()->global_min_z() > znew
      || uvw->domain()->global_max_z() < znew ) {
      boil::oout<<"pathline_advance:WARNING!!! Particle(no."<<ip
                <<") is outside of the domain. x,y,z= "
                <<xold<<" "<<yold<<" "<<zold<<"\n";
      particles[ip].u(0.0);
      particles[ip].v(0.0);
      particles[ip].w(0.0);
      xnew = xold;
      ynew = yold;
      znew = zold;
    }
    // store position
    particles[ip].x(xnew);
    particles[ip].y(ynew);
    particles[ip].z(znew);
  }

  boil::timer.stop("solidparticle advance");

  return;
}


