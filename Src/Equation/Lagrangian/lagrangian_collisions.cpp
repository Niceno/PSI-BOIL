#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::collisions() {
/*-------------------------------------------------+
|  check particle-particle and particle-wall collisions  |
+-------------------------------------------------*/

  boil::oout << "@lagrangian_collision.. " << boil::endl; 

  /*--------------------------------------------------------------------------+
  |  Browse through all the particle and calculate the minimum time for every  |
  |  particle-particle and particle-wall contact. Then check if this time is        |
  |  smaller than then time step ...                                          |
  |  if yes, advance the particle and take into consideration the collision    |
  +--------------------------------------------------------------------------*/
  const int nb_proc = boil::cart.nproc();

  boil::timer.start("lagrangian collisions");

  real nwall_x, nwall_y, nwall_z; /* normal vector components to the wall */
  real cur_time = 0.0;  /* current time: between zero and time step dt */

  while (cur_time < time->dt()) {
  
    loop_start_again: 

    cell_link();

    std::vector<int> indice_pa_wall; /* contains indices of particles
                                               colliding with wall */
    std::vector<double> wall_normal; /* contains which wall the
                                      particles are colliding with */
    std::vector<int> indice_p; /* contains indices of particle-particle 
                                                        collisions */
    real min_dt = time->dt(); /* min time of all particle-particle and 
                                 particle-wall collision */

    real min_dt_wall = time->dt(); /* min time of all particle
                                        colliding with wall */
    int col_wall     = OFF; /* 1 for particle-wall   collision */
    int col_particle = OFF; /* 1 for particle-particle collision */

     /* browse through all the particle */
    for(int pa=0; pa < size(); pa++) {

      /*------------------------------------------+
      |  1. First check particle-wall collisions  |
      +------------------------------------------*/
      /* check particle collision with wall */

      real dt_a_wall; /* time for particle a to collide with the wall */

      for (int b=0; b < u->bc(Comp::u()).count(); b++) {
        //if() {}
      }

      /*---------------------------------------------+
      |  2. Then check particle-particle collisions  |
      +---------------------------------------------*/
      const real xp = particles[pa].x();
      const real yp = particles[pa].y();
      const real zp = particles[pa].z();

      /* (ic,jc,kc) are the cell indices (in linked list) of pa */
      int ic = int (trunc ( ( xp - dom->global_min_x() ) / diam_x));
      int jc = int (trunc ( ( yp - dom->global_min_y() ) / diam_y));
      int kc = int (trunc ( ( zp - dom->global_min_z() ) / diam_z));

      if(ic == NX_coarse) ic--;
      if(jc == NY_coarse) jc--;
      if(kc == NZ_coarse) kc--;
      
      int west   = ic -1; int east  = ic +1;
      int south  = jc -1; int north = jc +1;
      int bottom = kc -1; int top   = kc +1;
 
      /* take care of boundary conditions */
      if(west   ==     OFF)  west    += 1;
      if(east   == NX_coarse)  east  -= 1;
      if(south  ==     OFF)  south   += 1;
      if(north  == NY_coarse)  north -= 1;
      if(bottom ==     OFF)  bottom  += 1;
      if(top    == NZ_coarse)  top   -= 1;

      /* browse through the neighboring cells */
      for (int ii = west; ii <= east; ii++) {
        //if() {}
      }  /* browsing in the close region */

    }  /* pa */

   delete [] link;


   if(col_wall == 1 && col_particle == 1) {
     if(min_dt > min_dt_wall) {
       min_dt = min_dt_wall;
       col_particle = 0;
     } else if (min_dt < min_dt_wall) {
       col_wall = 0; 
     }
   } else if (col_wall == 1) {
     min_dt = min_dt_wall;
   } 
   assert(min_dt <= time->dt());

   cur_time += min_dt;

   if (cur_time <= time->dt()) {
     for_p(p)
       for_m(m) {
         particles[p].xyz(m) += min_dt * particles[p].uvw(m);

         boil::cart.sum_real(& particles[p].xyz(m));
         /* averaging through all the processor is very important.
            In some cases, very small differences can build up and 
            cause the code to wait forever                       */ 
         particles[p].xyz(m) /= (real)nb_proc; 
       }
   } else {
     for_p(p)
       for_m(m) {
         particles[p].xyz(m) += (time->dt() - cur_time + min_dt) *
                                              particles[p].uvw(m);  //check-real-dt

         boil::cart.sum_real(& particles[p].xyz(m));
         particles[p].xyz(m) /= (real)nb_proc;
       }
   }

   if(cur_time <= time->dt()) {

     /* particle wall collision */           
     if(col_wall == 1) {
       //boil::oout <<"WALL_COLLISION.. " << wsize << boil::endl; 
     } 

     /* particle-particle collision-coalescence */
     if (col_particle == 1) {
       //if() {}
     } /* collisions and coalescence */ 

     indice_pa_wall.clear();
     wall_normal.clear();
     indice_p.clear();

   } /* cur_time <= dt */


   /*=========================================+
   |                                          |
   |  erase particle if it leaves the domain  |
   |                                          |
   +=========================================*/

   //loop_again: 

   for_p(p) { 
     const real xp = particles[p].x();
     const real yp = particles[p].y();
     const real zp = particles[p].z();
     const real r  = 0.5 * particles[p].d();

     int ind = 0;
      
     //for (int b=0; b < u->bc(Comp::u()).count(); b++) {
       //if(u->bc(Comp::u()).type(b) == BndType::outlet()) {

         //Dir d = u->bc(Comp::u()).direction(b);

         //if(d == Dir::imax() && fabs(xp - dom->global_max_x()) <= r) 
         //if(d == Dir::imax() && (xp >= dom->global_max_x()) ) 
         //if( xp >= dom->global_max_x() )
         if( xp >= (dom->global_max_x()-r) ) {
           OPR("Leaving imax.."); erase(p); ind = 1;
           std::ofstream outfile("erase-particle-after-leaving-xmax.txt", std::ios_base::app);  //appendix-write-in
           outfile << xp << boil::endl;
           outfile.close();
         }
         //if(d == Dir::imin() && fabs(xp - dom->global_min_x()) <= r) 
         if( xp <= (dom->global_min_x()+r) ) {
           OPR("Leaving imin.."); erase(p); ind = 1;
           std::ofstream outfile("erase-particle-after-leaving-xmin.txt", std::ios_base::app);  //appendix-write-in
           outfile << xp << boil::endl;
           outfile.close();
         }
         //if(d == Dir::jmax() && fabs(yp - dom->global_max_y()) <= r) 
         if( yp >= (dom->global_max_y()-r) ) {
           OPR("Leaving jmax.."); erase(p); ind = 1;
           std::ofstream outfile("erase-particle-after-leaving-ymax.txt", std::ios_base::app);  //appendix-write-in
           outfile << yp << boil::endl;
           outfile.close();
         }
         //if(d == Dir::jmin() && fabs(yp - dom->global_min_y()) <= r) 
         //if(d == Dir::jmin() && fabs(yp - dom->global_min_y()) <= r) 
         //if(d == Dir::jmin() && (yp <= dom->global_min_y()) ) 
         //if( yp <= dom->global_min_y() )
         if( yp <= (dom->global_min_y()+r) ) {
           OPR("Leaving jmin.."); erase(p); ind = 1;
           std::ofstream outfile("erase-particle-after-leaving-ymin.txt", std::ios_base::app);  //appendix-write-in
           outfile << yp << boil::endl;
           outfile.close();
         }
         //if(d == Dir::kmax() && fabs(zp - dom->global_max_z()) <= r) 
         //if( zp >= (dom->global_max_z()-r) )  //error-over-defined-with-color-function-criterion
         if( zp >= (dom->global_max_z()-r) ) {
           OPR("Leaving kmax.."); erase(p); ind = 1;
           std::ofstream outfile("erase-particle-after-leaving-zmax.txt", std::ios_base::app);  //appendix-write-in
           outfile << zp << boil::endl;
           outfile.close();
         }
         //if(d == Dir::kmin() && fabs(zp - dom->global_min_z()) <= r) 
         if( zp <= (dom->global_min_z()+r) ) {
           OPR("Leaving kmin.."); erase(p); ind = 1;
           std::ofstream outfile("erase-particle-after-leaving-zmin.txt", std::ios_base::app);  //appendix-write-in
           outfile << zp << boil::endl;
           outfile.close();
         }

         //check-color-function-to-"erase"/"deposit"-particles
         //define-the-interface-of-continuous-phases
         const int cell_i = ( dom->local_i( dom->I(particles[p].x()) ) );
         const int cell_j = ( dom->local_j( dom->J(particles[p].y()) ) );
         const int cell_k = ( dom->local_k( dom->K(particles[p].z()) ) );
         if( (*cfu)[cell_i][cell_j][cell_k] > 0.5 ) {
           OPR("Deposition.."); erase(p); ind = 1;
           std::ofstream outfile("deposit-particle-into-liquid.txt", std::ios_base::app);  //appendix-write-in
           outfile << xp << " " << yp << " " << zp << boil::endl;
           outfile.close();
         }

       //}
     //}
     //if (ind == 1) {goto loop_again;} 
   }
 }
 
 boil::timer.stop("lagrangian collisions");

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_collisions.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
