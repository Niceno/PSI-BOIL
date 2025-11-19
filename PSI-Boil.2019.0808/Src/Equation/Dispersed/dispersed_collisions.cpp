#include "dispersed.h"

/******************************************************************************/
void Dispersed::collisions() {
/*-------------------------------------------------+
|  check bubble-bubble and bubble-wall collisions  |
+-------------------------------------------------*/

  boil::oout << "@dispersed_collision.. " << boil::endl; 

  /*--------------------------------------------------------------------------+
  |  Browse through all the bubbles and calculate the minimum time for every  |
  |  bubble-bubble and bubble-wall contact. Then check if this time is        |
  |  smaller than then time step ...                                          |
  |  if yes, advance the bubbles and take into consideration the collision    |
  +--------------------------------------------------------------------------*/
  const real sigma    = flu->sigma()->value(); 
  const real mu_cont  = flu->mu(continuous); /* viscosity carrier phase */
  const real rho_disp = flu->rho(dispersed);
  const real rho_cont = flu->rho(continuous);

  const int nb_proc = boil::cart.nproc();

  boil::timer.start("dispersed collisions");

  real nwall_x, nwall_y, nwall_z; /* normal vector components to the wall */
  real cur_time = 0.0;  /* cuurent time: between zero and time step dt */

  while (cur_time < time->dt()) {
  
    loop_start_again: 

    cell_link();

    std::vector<int> indice_pa_wall; /* contains indices of particles
                                               colliding with wall */
    std::vector<double> wall_normal; /* contains which wall the
                                      particles are colliding with */
    std::vector<int> indice_p; /* contains indices of bubble-bublle 
                                                        collisions */
    real min_dt = time->dt(); /* min time of all bubble-bubble and 
                                 bubble-wall collision */

    real min_dt_wall = time->dt(); /* min time of all bubbles
                                        colliding with wall */
    int col_wall    = OFF; /* 1 for bubble-wall   collision */
    int col_bubbles = OFF; /* 1 for bubble-bubble collision */

     /* browse through all the bubbles */
    for(int pa=0; pa < size(); pa++) { 

      /* check bubble collision with wall */

      real dt_a_wall; /* time for bubble a to collide with the wall */

      for (int b=0; b < u->bc(Comp::u()).count(); b++) {

        if(u->bc(Comp::u()).type(b) == BndType::wall()) {

          Dir d = u->bc(Comp::u()).direction(b);
          
     	  /*-------+
          |  imin  |
          +-------*/
          if(d == Dir::imin()) {
            real dis_aw  = dom->global_min_x() - particles[pa].x();
            real ua = particles[pa].uvw(Comp::u());             
            if(ua < 0.0 && fabs(dis_aw) <= 3.0 * particles[pa].d()) {
              dt_a_wall = (0.5 * particles[pa].d() +  dis_aw) / ua;

              if(dt_a_wall <= min_dt_wall) {

                if(dt_a_wall < min_dt_wall) {
                  indice_pa_wall.clear();               
                  wall_normal.clear();
                }
                min_dt_wall = dt_a_wall;
                indice_pa_wall.push_back (pa);               
                wall_normal.push_back(1.0);
                wall_normal.push_back(0.0);
                wall_normal.push_back(0.0);
                col_wall  = 1;
              }
            }
          }
 
     	  /*-------+
          |  imax  |
          +-------*/
          if(d == Dir::imax()) {
            real dis_aw  = particles[pa].x() - dom->global_max_x();
            real ua = -particles[pa].uvw(Comp::u());             
            if(ua < 0.0 && fabs(dis_aw) <= 3.0 * particles[pa].d()) {
              dt_a_wall = (0.5 * particles[pa].d() +  dis_aw) / ua;

              if(dt_a_wall <= min_dt_wall) {
                if(dt_a_wall < min_dt_wall) {
                  wall_normal.clear();
                  indice_pa_wall.clear();               
                }
                min_dt_wall = dt_a_wall;
                indice_pa_wall.push_back (pa);               
                wall_normal.push_back(-1.0);
                wall_normal.push_back(0.0);
                wall_normal.push_back(0.0);
                col_wall  = 1;
              }
            }
          }

     	  /*-------+
          |  jmin  |
          +-------*/
          if(d == Dir::jmin()) {
            real dis_aw  = dom->global_min_y() - particles[pa].y();
            real ua = particles[pa].uvw(Comp::v());             
            if(ua < 0.0 && fabs(dis_aw) <= 3.0 * particles[pa].d()) {
              dt_a_wall = (0.5 * particles[pa].d() +  dis_aw) / ua;

              if(dt_a_wall <= min_dt_wall) {
                if(dt_a_wall < min_dt_wall) {
                  indice_pa_wall.clear();               
                  wall_normal.clear();
                }
                min_dt_wall = dt_a_wall;
                indice_pa_wall.push_back (pa);               
                wall_normal.push_back(0.0);
                wall_normal.push_back(1.0);
                wall_normal.push_back(0.0);
                col_wall  = 1;

              }
            }
          }
         
     	  /*-------+
          |  jmax  |
          +-------*/
          if(d == Dir::jmax()) {
            real dis_aw  = particles[pa].y() - dom->global_max_y();
            real ua = -particles[pa].uvw(Comp::v());             
            if(ua < 0.0 && fabs(dis_aw) <= 3.0 * particles[pa].d()) {
              dt_a_wall = (0.5 * particles[pa].d() +  dis_aw) / ua;
              if(dt_a_wall <= min_dt_wall) {
                if(dt_a_wall < min_dt_wall) {
                  indice_pa_wall.clear();               
                  wall_normal.clear();
                }
                min_dt_wall = dt_a_wall;
                indice_pa_wall.push_back (pa);               
                wall_normal.push_back(0.0);
                wall_normal.push_back(-1.0);
                wall_normal.push_back(0.0);
                col_wall  = 1;
              }
            }
          }

     	  /*-------+
          |  kmin  |
          +-------*/
          if(d == Dir::kmin()) {
            real dis_aw  = dom->global_min_z() - particles[pa].z();
            real ua = particles[pa].uvw(Comp::w());             
            if(ua < 0.0 && fabs(dis_aw) <= 3.0 * particles[pa].d()) {
              dt_a_wall = (0.5 * particles[pa].d() +  dis_aw) / ua;

              if(dt_a_wall <= min_dt_wall) {
                if(dt_a_wall < min_dt_wall) {
                 indice_pa_wall.clear();               
                 wall_normal.clear();
                }
                min_dt_wall = dt_a_wall;
                indice_pa_wall.push_back (pa);               
                wall_normal.push_back(0.0);
                wall_normal.push_back(0.0);
                wall_normal.push_back(1.0);
                col_wall  = 1;
              }
            }
          }

     	  /*-------+
          |  kmax  |
          +-------*/
          if(d == Dir::kmax()) {
            real dis_aw  = particles[pa].z() - dom->global_max_z();
            real ua = -particles[pa].uvw(Comp::w());             
            if(ua < 0.0 && fabs(dis_aw) <= 3.0 * particles[pa].d()) {
              dt_a_wall = (0.5 * particles[pa].d() +  dis_aw) / ua;

              if(dt_a_wall <= min_dt_wall) {
                if(dt_a_wall < min_dt_wall) {
                  indice_pa_wall.clear();               
                  wall_normal.clear();
                }
                min_dt_wall = dt_a_wall;
                indice_pa_wall.push_back (pa);               
                wall_normal.push_back(0.0);
                wall_normal.push_back(0.0);
                wall_normal.push_back(-1.0);
                col_wall  = 1;
              }
            }
          }
        }
      }
     
      /*-------------------------------------+
      |  now check bubble-bubble collisions  |
      +-------------------------------------*/
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
        for (int jj = south; jj <= north; jj++) {
          for (int kk = bottom; kk <= top; kk++) {

            /* start from the head of the chain */
            int pb = cell[ii][jj][kk];
            
            while (pb > pa) {
              const real dx = particles[pa].x() - particles[pb].x();
              const real dy = particles[pa].y() - particles[pb].y();
              const real dz = particles[pa].z() - particles[pb].z();
              const real dist_ab  = sqrt(dx*dx + dy*dy + dz*dz);
              const real radius_a = particles[pa].d() * 0.5;
              const real radius_b = particles[pb].d() * 0.5;

              real rel_dis = dist_ab - (radius_a + radius_b - boil::pico);
              if(rel_dis < 0.0) OPR(rel_dis);
 
              real r_v = dx * (particles[pa].uvw(Comp::u())  - 
                               particles[pb].uvw(Comp::u())) +
                         dy * (particles[pa].uvw(Comp::v())  - 
                               particles[pb].uvw(Comp::v())) +
                         dz * (particles[pa].uvw(Comp::w())  - 
                               particles[pb].uvw(Comp::w())) ;

              /* check if particles intersect eachother */
              if ( dist_ab < boil::maxr(radius_a,radius_b) ) {
                boil::oout <<"bubbles inside eachother..merge them.. " << pa << 
                                                       " " << pb << boil::endl;
                merge(pa,pb);
                goto loop_start_again;
              }

              const real r_v_2 = r_v * r_v;
                    real v_2   = (particles[pa].uvw(Comp::u())  -
                                  particles[pb].uvw(Comp::u())) *  
                                 (particles[pa].uvw(Comp::u())  -
                                  particles[pb].uvw(Comp::u())) + 

                                 (particles[pa].uvw(Comp::v())  -
                                  particles[pb].uvw(Comp::v())) *  
                                 (particles[pa].uvw(Comp::v())  -
                                  particles[pb].uvw(Comp::v())) + 
 
                                 (particles[pa].uvw(Comp::w())  -
                                  particles[pb].uvw(Comp::w())) *  
                                 (particles[pa].uvw(Comp::w())  -
                                  particles[pb].uvw(Comp::w())) ;

              const real r_2  = dist_ab * dist_ab;
              const real radi =  (radius_a + radius_b); 
              real discriminant = r_v_2 - v_2 *(r_2 - radi * radi); 
 
              /* bubbles approaching each other */
              if (r_v < 0.0 && discriminant >= 0.0 && v_2 > 0.0) {

                real dt_ab =(-r_v - sqrt(discriminant)) / v_2; 

                if(dt_ab < 0.0)  dt_ab = 0.0;

                assert(dt_ab >= 0.0);
 
                if(dt_ab <= min_dt) {
                  if(dt_ab < min_dt) indice_p.clear();               
                  min_dt = dt_ab;
                  col_bubbles = 1;
                  indice_p.push_back (pa);
                  indice_p.push_back (pb);
                }                     
 
              } /* moving towards each other */

              pb = link[pb];

            }  /* while pb != OFF */
          }
        }
      }  /* browsing in the close region */
    }  /* pa */

   delete [] link;

   if(col_wall == 1 && col_bubbles == 1) {
     if(min_dt > min_dt_wall) {
       min_dt = min_dt_wall;
       col_bubbles = 0;
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
                                              particles[p].uvw(m);
          boil::cart.sum_real(& particles[p].xyz(m));
          particles[p].xyz(m) /= (real)nb_proc;
       }
   }

   if(cur_time <= time->dt()) {

     /* bubble wall collision */           
     if(col_wall == 1) { 
       int wsize = indice_pa_wall.size();
       boil::oout <<"WALL_COLLISION.. " << wsize << boil::endl; 

       for(int a_wall = 0; a_wall < indice_pa_wall.size(); a_wall++) { 
         int pa_wall = indice_pa_wall[a_wall];
         real ux = particles[pa_wall].uvw(Comp::u());
         real uy = particles[pa_wall].uvw(Comp::v());
         real uz = particles[pa_wall].uvw(Comp::w());
         real nwall_x = wall_normal[3*a_wall];
         real nwall_y = wall_normal[3*a_wall +1];
         real nwall_z = wall_normal[3*a_wall +2];

         real un = fabs(ux * nwall_x + uy * nwall_y + uz * nwall_z);
         particles[pa_wall].uvw(Comp::u()) += 2.0 * un * nwall_x;
         particles[pa_wall].uvw(Comp::v()) += 2.0 * un * nwall_y;
         particles[pa_wall].uvw(Comp::w()) += 2.0 * un * nwall_z;
       }
     } 

     /* bubble-bubble collision-coalescence */
     if (col_bubbles == 1) {
       int psize = indice_p.size()/2;
       for(int ab = 0; ab < indice_p.size()/2; ab++) { 
         int indice_a = indice_p[2*ab];
         int indice_b = indice_p[2*ab+1];
         real dx = particles[indice_a].x() - particles[indice_b].x();
         real dy = particles[indice_a].y() - particles[indice_b].y();
         real dz = particles[indice_a].z() - particles[indice_b].z();
         real dist_ab  = sqrt(dx*dx + dy*dy + dz*dz);
         dx /= dist_ab; dy /= dist_ab; dz /= dist_ab;
          
         
         #if COALESCENCE == false   /* consider only rebouncing */
         rebounce(indice_a,indice_b); 
         #else /* consider coalescence and rebouncing */
         real ua_n = particles[indice_a].uvw(Comp::u()) * dx +   
                     particles[indice_a].uvw(Comp::v()) * dy +   
                     particles[indice_a].uvw(Comp::w()) * dz ;   
         real ub_n = particles[indice_b].uvw(Comp::u()) * dx +   
                     particles[indice_b].uvw(Comp::v()) * dy +   
                     particles[indice_b].uvw(Comp::w()) * dz ;   

         real r_eqv = 1.0 / (1.0 / particles[indice_a].d() +
                             1.0 / particles[indice_b].d() );    
         real u_rel = fabs(ua_n - ub_n);
         assert (u_rel > 0.0);

         /* contact time calculated from Sommerfeld model.
            0.5 is the calibration factor */
          real contact_time  = 0.5  * r_eqv / u_rel;

         /* drainage time calculated from Prince and Blanch model 
            9.21 = ln(h0/hf) with h0 = 10^-4 and hf = 10^-8 */
         real drainage_time = 9.21 * sqrt(r_eqv * r_eqv * r_eqv *
                                       rho_cont /(16.0 * sigma));
 
         if(drainage_time > contact_time) { 
           OPR("BUBBLE-BUBBLE REBOUNCING ");
           rebounce(indice_a,indice_b); 
         } else {
           OPR("BUBBLE-BUBBLE COALESCENCE");
           merge(indice_a,indice_b); 
           goto loop_start_again;
         } 
         #endif
       }
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

   loop_again: 

   for_p(p) { 
     const real xp = particles[p].x();
     const real yp = particles[p].y();
     const real zp = particles[p].z();
     const real r  = 0.5 * particles[p].d();

     int ind = 0;
      
     for (int b=0; b < u->bc(Comp::u()).count(); b++) {
       if(u->bc(Comp::u()).type(b) == BndType::outlet()) {

         Dir d = u->bc(Comp::u()).direction(b);

         if(d == Dir::imax() && fabs(xp - dom->global_max_x()) <= r) 
           {OPR("Leaving imax.."); erase(p);ind = 1;}
         if(d == Dir::imin() && fabs(xp - dom->global_min_x()) <= r) 
           {OPR("Leaving imin.."); erase(p);ind = 1;}
         if(d == Dir::jmax() && fabs(yp - dom->global_max_y()) <= r) 
           {OPR("Leaving jmax.."); erase(p);ind = 1;}
         if(d == Dir::jmin() && fabs(yp - dom->global_min_y()) <= r) 
           {OPR("Leaving jmin.."); erase(p);ind = 1;}
         if(d == Dir::kmax() && fabs(zp - dom->global_max_z()) <= r) 
           {OPR("Leaving kmax.."); erase(p);ind = 1;}
         if(d == Dir::kmin() && fabs(zp - dom->global_min_z()) <= r) 
           {OPR("Leaving kmin.."); erase(p);ind = 1;}
       }
     }
     if (ind == 1) {goto loop_again;} 
   }
 }
 
 boil::timer.stop("dispersed collisions");

}
