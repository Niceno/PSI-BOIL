#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::forces() {
/*=================================================+
|                                                  |
|  browse through all particles to compute forces  |
|                                                  |
+=================================================*/

  const int nb_proc = boil::cart.nproc();

  /* morton number */
  /*const real morton = boil::g * delta_rho *
                      mu_cont * mu_cont * mu_cont * mu_cont /
                      (rho_cont * rho_cont * sigma * sigma * sigma);
  const real coef_morton = boil::maxr(6.0 * log10(morton) + 24.0, 4.4);*/  //not-in-use

  //gravity-?!


  for_p(p) { /* browse through all particles */
      
    /*--------------------------------------+
    |  browse through all corners of a box  |
    |                                       |
    |            6---------------7          |
    |           /|      T       /|          |
    |       ke 4---------------5 |          |
    |          | |       N     | |          |
    |          |W|      /      |E|          |
    |          | |     S       | |          |
    |          | 2-------------|-3 je       |
    |          |/       r_y      |/           |
    |       ks 0---------------1 js         |
    |         is              ie            |
    +--------------------------------------*/
	
    /* velocities at box's faces */
    real uvw_w[DIM] = {0.0, 0.0, 0.0};
    real uvw_e[DIM] = {0.0, 0.0, 0.0};
    real uvw_s[DIM] = {0.0, 0.0, 0.0};
    real uvw_n[DIM] = {0.0, 0.0, 0.0};
    real uvw_b[DIM] = {0.0, 0.0, 0.0};
    real uvw_t[DIM] = {0.0, 0.0, 0.0};

    for_m(m) {
      uvw_w[~m] = 0.25 * (particles[p].box_uvw(0,m)+particles[p].box_uvw(2,m) +
                          particles[p].box_uvw(4,m)+particles[p].box_uvw(6,m));  
      uvw_e[~m] = 0.25 * (particles[p].box_uvw(1,m)+particles[p].box_uvw(3,m) +
                          particles[p].box_uvw(5,m)+particles[p].box_uvw(7,m)); 
      uvw_s[~m] = 0.25 * (particles[p].box_uvw(0,m)+particles[p].box_uvw(1,m) +
                          particles[p].box_uvw(4,m)+particles[p].box_uvw(5,m));  
      uvw_n[~m] = 0.25 * (particles[p].box_uvw(2,m)+particles[p].box_uvw(3,m) +
                          particles[p].box_uvw(6,m)+particles[p].box_uvw(7,m));  
      uvw_b[~m] = 0.25 * (particles[p].box_uvw(0,m)+particles[p].box_uvw(1,m) +
                          particles[p].box_uvw(2,m)+particles[p].box_uvw(3,m));  
      uvw_t[~m] = 0.25 * (particles[p].box_uvw(4,m)+particles[p].box_uvw(5,m) +
                          particles[p].box_uvw(6,m)+particles[p].box_uvw(7,m));  
    } 

    int BoxingScheme = 2;  //interpolation-mark-modify-into-main.cpp
    real uvwc[DIM] = {0.0, 0.0, 0.0};

    real mu_cont   = -99.1;  //can't be const, and set initial value for checking
    real rho_cont  = -99.2;
    real rho_lagp  = -99.3;
    real delta_rho = -99.4;

    if (BoxingScheme==1) {
      //step-1-locate particle position and achieve cell coordinate
      const int cell_i = ( dom->local_i( dom->I(particles[p].x()) ) );
      const int cell_j = ( dom->local_j( dom->J(particles[p].y()) ) );
      const int cell_k = ( dom->local_k( dom->K(particles[p].z()) ) );

      //calculate-the-properties-of-all-the-three-phases
      //continuous
      //lagparticle
      mu_cont  =  flu->mu(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                + flu->mu(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //1-mu_cont
      rho_cont =  flu->rho(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                           + flu->rho(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //2-rho_cont
      rho_lagp = sol->rho(cell_i,cell_j,cell_k);  //3-rho_lagp-const
      delta_rho = fabs(rho_cont  - rho_lagp);  //4-delta_rho

      /*const real mu_cont  =  flu->mu(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )  //error
                           + flu->mu(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //1-mu_cont
      const real rho_cont =  flu->rho(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                           + flu->rho(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //2-rho_cont
      const real rho_lagp = sol->rho(cell_i,cell_j,cell_k);  //3-rho_lagp-const
      const real delta_rho = fabs(rho_cont  - rho_lagp);  //4-delta_rho */

      /* undisturbed liquid velocity */
      /*uvwc[DIM] = {0.5 * (uvw_w[0] + uvw_e[0]), 
                   0.5 * (uvw_s[1] + uvw_n[1]), 
                   0.5 * (uvw_b[2] + uvw_t[2])}; }*/
      uvwc[0] = 0.5 * (uvw_w[0] + uvw_e[0]);
      uvwc[1] = 0.5 * (uvw_s[1] + uvw_n[1]);
      uvwc[2] = 0.5 * (uvw_b[2] + uvw_t[2]);
      }
    else if (BoxingScheme==2) {
      //step-1-locate particle position and achieve cell coordinate
      const int cell_i = ( dom->local_i( dom->I(particles[p].x()) ) );
      const int cell_j = ( dom->local_j( dom->J(particles[p].y()) ) );
      const int cell_k = ( dom->local_k( dom->K(particles[p].z()) ) );

      //calculate-the-properties-of-all-the-three-phases
      //continuous
      //lagp
      mu_cont  =  flu->mu(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                + flu->mu(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //1-mu_cont
      rho_cont =  flu->rho(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                + flu->rho(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //2-rho_cont
      rho_lagp = sol->rho(cell_i,cell_j,cell_k);  //3-rho_lagp-const
      delta_rho = fabs(rho_cont  - rho_lagp);  //4-delta_rho

      //step-2.2-interpolate continuous variables on particle position
      uvwc[0] = (*u)[Comp::u()][cell_i][cell_j][cell_k];
      uvwc[1] = (*u)[Comp::v()][cell_i][cell_j][cell_k];
      uvwc[2] = (*u)[Comp::w()][cell_i][cell_j][cell_k];
      }
    else if (BoxingScheme==3) {
      //step-1-locate particle position and achieve cell coordinate
      const int cell_i = ( dom->local_i( dom->I(particles[p].x()) ) );
      const int cell_j = ( dom->local_j( dom->J(particles[p].y()) ) );
      const int cell_k = ( dom->local_k( dom->K(particles[p].z()) ) );

      //calculate-the-properties-of-all-the-three-phases
      //continuous
      //lagp
      mu_cont  =  flu->mu(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                + flu->mu(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //1-mu_cont
      rho_cont =  flu->rho(1.0) * ( (*cfu)[cell_i][cell_j][cell_k] )
                + flu->rho(0.0) * ( 1.0 - (*cfu)[cell_i][cell_j][cell_k] );  //2-rho_cont
      rho_lagp = sol->rho(cell_i,cell_j,cell_k);  //3-rho_lagp-const
      delta_rho = fabs(rho_cont  - rho_lagp);  //4-delta_rho

      //step-2.2-interpolate continuous variables on particle position
      /*if( cell_i==1 || cell_i==(ni()-2)
          || cell_j==1 || cell_j==(nj()-2)
          || cell_k==1 || cell_k==(nk()-2) ) {
        uvwc[0] = (*u)[Comp::u()][cell_i][cell_j][cell_k];
        uvwc[1] = (*u)[Comp::v()][cell_i][cell_j][cell_k];
        uvwc[2] = (*u)[Comp::w()][cell_i][cell_j][cell_k]; } */
      if( cell_j==1 ) {
        //cells-at-wall
        //uvwc[0] = (*u)[Comp::u()][cell_i][cell_j][cell_k]/(2.5e-5)*particles[p].y();
        //uvwc[1] = (*u)[Comp::v()][cell_i][cell_j][cell_k]/(2.5e-5)*particles[p].y();
        uvwc[0] = (*u)[Comp::u()][cell_i][cell_j][cell_k];
        uvwc[1] = (*u)[Comp::v()][cell_i][cell_j][cell_k];
        uvwc[2] = (*u)[Comp::w()][cell_i][cell_j][cell_k];
        }
      else {
        //3D-independent-2nd-order-Lagrange-interpolation     
        real LX[3],LY[3],LZ[3];

        //step-2.3.u-interpolation
        LX[0] = (particles[p].x() - xc(cell_i)) / (xc(cell_i-1) - xc(cell_i))
               *(particles[p].x() - xc(cell_i+1)) / (xc(cell_i-1) - xc(cell_i+1));
        LX[1] = (particles[p].x() - xc(cell_i-1)) / (xc(cell_i) - xc(cell_i-1))
               *(particles[p].x() - xc(cell_i+1)) / (xc(cell_i) - xc(cell_i+1));
        LX[2] = (particles[p].x() - xc(cell_i-1)) / (xc(cell_i+1) - xc(cell_i-1))
               *(particles[p].x() - xc(cell_i)) / (xc(cell_i+1) - xc(cell_i));
        uvwc[0] = (LX[0] * 1.0 * 1.0) * (*u)[Comp::u()][cell_i-1][cell_j][cell_k]
                 +(LX[1] * 1.0 * 1.0) * (*u)[Comp::u()][cell_i][cell_j][cell_k]
                 +(LX[2] * 1.0 * 1.0) * (*u)[Comp::u()][cell_i+1][cell_j][cell_k];

        //step-2.3.v-interpolation
        LY[0] = (particles[p].y() - yc(cell_j)) / (yc(cell_j-1) - yc(cell_j))
               *(particles[p].y() - yc(cell_j+1)) / (yc(cell_j-1) - yc(cell_j+1));
        LY[1] = (particles[p].y() - yc(cell_j-1)) / (yc(cell_j) - yc(cell_j-1))
               *(particles[p].y() - yc(cell_j+1)) / (yc(cell_j) - yc(cell_j+1));
        LY[2] = (particles[p].y() - yc(cell_j-1)) / (yc(cell_j+1) - yc(cell_j-1))
               *(particles[p].y() - yc(cell_j)) / (yc(cell_j+1) - yc(cell_j));
        uvwc[1] = (1.0 * LY[0] * 1.0) * (*u)[Comp::v()][cell_i][cell_j-1][cell_k]
                 +(1.0 * LY[1] * 1.0) * (*u)[Comp::v()][cell_i][cell_j][cell_k]
                 +(1.0 * LY[2] * 1.0) * (*u)[Comp::v()][cell_i][cell_j+1][cell_k];

        //step-2.3.w-interpolation
        LZ[0] = (particles[p].z() - zc(cell_k)) / (zc(cell_k-1) - zc(cell_k))
               *(particles[p].z() - zc(cell_k+1)) / (zc(cell_k-1) - zc(cell_k+1));
        LZ[1] = (particles[p].z() - zc(cell_k-1)) / (zc(cell_k) - zc(cell_k-1))
               *(particles[p].z() - zc(cell_k+1)) / (zc(cell_k) - zc(cell_k+1));
        LZ[2] = (particles[p].z() - zc(cell_k-1)) / (zc(cell_k+1) - zc(cell_k-1))
               *(particles[p].z() - zc(cell_k)) / (zc(cell_k+1) - zc(cell_k));
        uvwc[2] = (1.0 * 1.0 * LZ[0]) * (*u)[Comp::w()][cell_i][cell_j][cell_k-1]
                 +(1.0 * 1.0 * LZ[1]) * (*u)[Comp::w()][cell_i][cell_j][cell_k]
                 +(1.0 * 1.0 * LZ[2]) * (*u)[Comp::w()][cell_i][cell_j][cell_k+1];
        }
    }


    /*real dudx = (uvw_e[0] - uvw_w[0]) / particles[p].box_dx(); 
    real dudy = (uvw_n[0] - uvw_s[0]) / particles[p].box_dy(); 
    real dudz = (uvw_t[0] - uvw_b[0]) / particles[p].box_dz(); 
    real dvdx = (uvw_e[1] - uvw_w[1]) / particles[p].box_dx(); 
    real dvdy = (uvw_n[1] - uvw_s[1]) / particles[p].box_dy(); 
    real dvdz = (uvw_t[1] - uvw_b[1]) / particles[p].box_dz(); 
    real dwdx = (uvw_e[2] - uvw_w[2]) / particles[p].box_dx(); 
    real dwdy = (uvw_n[2] - uvw_s[2]) / particles[p].box_dy(); 
    real dwdz = (uvw_t[2] - uvw_b[2]) / particles[p].box_dz();*/ 

    /*-----------------------------------------------+
    |  compute rot(u_liquid)= r_x i + r_y j + r_z k  |
    +------------------------------------------------+
    |         u_liquid = (u,v,w)                     |
    |         r_x = dw/dy - dv/dz                    |
    |         r_y = du/dz - dw/dx                    |
    |         r_z = dv/dx - du/dy                    |
    +-----------------------------------------------*/
    /*const real delta_x = particles[p].box_dx();
    const real delta_y = particles[p].box_dy();
    const real delta_z = particles[p].box_dz();
    const real r_x = (uvw_n[~Comp::w()] - uvw_s[~Comp::w()]) / delta_y 
                   - (uvw_t[~Comp::v()] - uvw_b[~Comp::v()]) / delta_z;   

    const real r_y = (uvw_t[~Comp::u()] - uvw_b[~Comp::u()]) / delta_z 
                   - (uvw_e[~Comp::w()] - uvw_w[~Comp::w()]) / delta_x;   

    const real r_z = (uvw_e[~Comp::v()] - uvw_w[~Comp::v()]) / delta_x 
                   - (uvw_n[~Comp::u()] - uvw_s[~Comp::u()]) / delta_y;*/   

    /* relative velocity */
    real uvwr[DIM] = {0.0, 0.0, 0.0}; 
    real norm_uvwr = 0.0;
    for_m(m) {
      uvwr[~m] = uvwc[~m] - (particles[p].uvw(m));
      norm_uvwr += uvwr[~m] * uvwr[~m];
    }
    norm_uvwr = sqrt(norm_uvwr);


    /*---------------------------------+
    |                                  |
    |  compute forces on the particle  |
    |                                  | 
    +---------------------------------*/
    real sum = 0.0; 
    real vol = 0.0;
    
    /* reynolds number based on relative velocity */
    const real reynolds = rho_cont * norm_uvwr * particles[p].d() 
                        / mu_cont;

    /* etvos number */
    /*const real eotvos = boil::g * delta_rho  
                      * particles[p].d() * particles[p].d() 
                      / sigma; */

    /*-----------+
    |  buoyancy  |
    +-----------*/
    /*real fb[DIM] = {0.0, 0.0, 0.0};*/

    /* assume bouyancy acts in positive z direction */
    /*fb[~Comp::w()] = (rho_cont - rho_lagp) * boil::g * particles[p].volume();*/
 
    /*-------------+
    |  added mass  |
    +-------------*/
    /*real fam[DIM] = {0.0, 0.0, 0.0};
    real cam = 0.5; /* added mass coefficient 
    real coef = rho_cont * (cam) * particles[p].volume();  // 1+cam //

    real uvw_gradu = uvwc[0] * dudx + uvwc[1] * dudy + uvwc[2] * dudz;
    real uvw_gradv = uvwc[0] * dvdx + uvwc[1] * dvdy + uvwc[2] * dvdz;
    real uvw_gradw = uvwc[0] * dwdx + uvwc[1] * dwdy + uvwc[2] * dwdz;

    real dudt = (uvwc[0] - particles[p].uvwc_old(Comp::u())) / time->dt();
    real dvdt = (uvwc[1] - particles[p].uvwc_old(Comp::v())) / time->dt();
    real dwdt = (uvwc[2] - particles[p].uvwc_old(Comp::w())) / time->dt();

    /* Du/Dt 
    real duvw_dt[DIM]= {dudt+uvw_gradu, dvdt+uvw_gradv,dwdt+uvw_gradw};

    particles[p].uvwc_old(Comp::u(), uvwc[0]);
    particles[p].uvwc_old(Comp::v(), uvwc[1]);
    particles[p].uvwc_old(Comp::w(), uvwc[2]);

    if(time->current_step() != 1) {    
      fam[~Comp::u()] = coef * duvw_dt[0];
      fam[~Comp::v()] = coef * duvw_dt[1];
      fam[~Comp::w()] = coef * duvw_dt[2];
    }*/
  
    /*-------+
    |  drag  |
    +-------*/
    real fd[DIM]   = {0.0, 0.0, 0.0};
    real cd = 0.44; /* drag coefficient */

    //#if TOMIYAMA_DRAG
    /*if (reynolds > 0.0) {
      const real a = 48.0 / reynolds;
      const real b = (16.0 / reynolds) * (1.0 + 0.15 * pow(reynolds,0.687));
      const real cd_pure_system = boil::minr(a,b);
      const real c = (8.0/3.0) * eotvos / (4.0 + eotvos);      
      cd = boil::maxr(cd_pure_system, c);
    } else {
      cd = 0.0;
      OPR("Warning: Reynolds number is 0");
    } 
    #endif */

    //particle-motion-drag-force
    if (reynolds <= 1000.0) {
      cd = 24.0/reynolds*(1.0+0.15*pow(reynolds,0.687));
    } else {
      cd = 0.44;
    }

    for_m(m) 
      fd[~m] = 0.5 * cd * rho_cont * particles[p].area() * norm_uvwr * uvwr[~m];

     
    /*-------+
    |  lift  |
    +-------*/
    /*real fl[DIM] = {0.0, 0.0, 0.0};
    real cl = 0.15; /* lift coefficient 

    #if TOMIYAMA_LIFT  
    real et_pow = pow(eotvos, 0.757);
    const real d_h = particles[p].d() * pow(1.0 + 0.163 * et_pow, 1.0/3.0);

    const real et_h = boil::g * (delta_rho) * d_h * d_h / sigma;
 
    const real f_eo = 0.00105 * et_h * et_h * et_h 
                    - 0.0159  * et_h * et_h
                    - 0.0204  * et_h 
                    + 0.474;
    real tan = 0.288 * tanh(0.121 * reynolds);

    if(et_h < 4) {
      cl = boil::minr(tan, f_eo);
    } else if (et_h <= 10) {
      cl = f_eo;
    } else {
      cl = -0.27;
    }
    #endif

    const real alpha_lift = cl * rho_cont * particles[p].volume();

    fl[0] = alpha_lift * (uvwr[1] * r_z - uvwr[2] * r_y); 
    fl[1] = alpha_lift * (uvwr[2] * r_x - uvwr[0] * r_z);  
    fl[2] = alpha_lift * (uvwr[0] * r_y - uvwr[1] * r_x);*/

    /*-------------------------+
    |  Wall lubrication force  |
    +-------------------------*/
    /*
    real fw[DIM] = {0.0, 0.0, 0.0};

    real channel_length_x  = dom->global_max_x() - dom->global_min_x();
    real channel_length_y  = dom->global_max_y() - dom->global_min_y();
    real channel_length_z  = dom->global_max_z() - dom->global_min_z();

    real cw; //* wall coef based on Hosokawa  model 2003 

    for (int b=0; b < u->bc(Comp::u()).count(); b++) {

      if(u->bc(Comp::u()).type(b) == BndType::wall()) {

        Dir d = u->bc(Comp::u()).direction(b);

        real db = particles[p].d(); 
        real dis_bw; 
         
        if(reynolds < boil::milli) {
          OPR("Warning: Re = 0.0.. set cw to zero ... ");
          cw = 0.0;
        } else {
          real cw1 = coef_morton / pow(reynolds, 1.9) ;
          cw = boil::maxr(cw1, 0.0217 * eotvos);
        }

        if(d == Dir::imin()) {
          dis_bw = fabs(particles[p].x() - dom->global_min_x());
          real dis_bw2 = channel_length_x - dis_bw;
          real ur_x = uvwr[1]*uvwr[1] + uvwr[2]*uvwr[2];             
          fw[0] += cw * 0.5 * particles[p].d() * rho_cont *
                   particles[p].volume() * ur_x * 
          ((1.0 / (dis_bw * dis_bw)) - (1.0 / (dis_bw2 * dis_bw2))); 
        }

        if(d == Dir::imax()) {
          real dis_bw  = fabs(particles[p].x() - dom->global_max_x());
          real dis_bw2 = channel_length_x - dis_bw;
          real ur_x = uvwr[1]*uvwr[1] + uvwr[2]*uvwr[2];
          fw[0] += -cw * 0.5 * particles[p].d() * rho_cont * 
                    particles[p].volume() * ur_x * 
          ((1.0 / (dis_bw * dis_bw)) - (1.0 / (dis_bw2 * dis_bw2))); 
        }

        if(d == Dir::jmin()) {
          dis_bw = fabs(particles[p].y() - dom->global_min_y());
          real dis_bw2 = channel_length_y - dis_bw;
          real ur_y = uvwr[0]*uvwr[0] + uvwr[2]*uvwr[2];
          fw[1] += cw * rho_cont * 0.5 * particles[p].d() * 
                   particles[p].volume() * ur_y * 
          ((1.0 / (dis_bw * dis_bw)) - (1.0 / (dis_bw2 * dis_bw2))); 
        }

        if(d == Dir::jmax()) {
          real dis_bw  = fabs(particles[p].y() - dom->global_max_y());
          real dis_bw2 = channel_length_y - dis_bw;
          real ur_y = uvwr[0]*uvwr[0] + uvwr[2]*uvwr[2];
          fw[1] += -cw * rho_cont * 0.5 * particles[p].d() * 
                    particles[p].volume() * ur_y * 
          ((1.0 / (dis_bw * dis_bw)) - (1.0 / (dis_bw2 * dis_bw2))); 
        }

        if(d == Dir::kmin()) {
          dis_bw = fabs(particles[p].z() - dom->global_min_z());
          real dis_bw2 = channel_length_z - dis_bw;
          real ur_z = uvwr[0]*uvwr[0] + uvwr[1]*uvwr[1];
          fw[2] += cw * rho_cont * 0.5 * particles[p].d() *  
                   particles[p].volume() * ur_z * 
          ((1.0 / (dis_bw * dis_bw)) - (1.0 / (dis_bw2 * dis_bw2))); 
        }

        if(d == Dir::kmax()) {
          real dis_bw  = fabs(particles[p].z() - dom->global_max_z());
          real dis_bw2 = channel_length_z - dis_bw;
          real ur_z = uvwr[0]*uvwr[0] + uvwr[1]*uvwr[1];
          fw[2] += -cw * rho_cont * 0.5 * particles[p].d() * 
                   particles[p].volume() * ur_z * 
          ((1.0 / (dis_bw * dis_bw)) - (1.0 / (dis_bw2 * dis_bw2))); 
        }
      }
    }
    */

    /*-------------+
    |  total force |
    +-------------*/
    real ft[DIM] = {0.0, 0.0, 0.0};

    for_m(m) 
      //ft[~m] = fb[~m] + fd[~m] + fl[~m] + fw[~m] + fam[~m]; 
      ft[~m] = fd[~m];  //only-drag-force-assumption

    /*---------------------------------------+
    |                                        |
    |  compute new velocity of the particle  |
    |                                        | 
    +---------------------------------------*/
    real duvw[DIM] = {0.0, 0.0, 0.0}; /* increase in velocity */

    for_m(m) {
      duvw[~m] = time->dt() * ft[~m] / (particles[p].volume() * (rho_lagp));
               //time->dt() * ft[~m] / (particles[p].volume() * (0.5 * rho_cont + rho_lagp));

      particles[p].uvw(m) += duvw[~m];

      /* averaging through all the processor is very important.
         In some cases, very small differences can build up and 
         cause the code to wait forever                       */ 
      boil::cart.sum_real(& particles[p].uvw(m));
      particles[p].uvw(m) /= (real)nb_proc;

      //if(particles[p].uvw(m) > uvw_limit) particles[p].uvw(m) = uvw_limit;
    }

  } /* browse through particles */

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_forces.cpp,v 1.2 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
