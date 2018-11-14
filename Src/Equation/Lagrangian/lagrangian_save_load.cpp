#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::save(const char * nm, const int it) {

  /* file name*/
  //std::string name = name_file(nm, ".bck", it);
  //std::string name = name_file(nm, ".dat", it);

  /* open a file */
  //std::ofstream out(name.c_str(), std::ios::binary);

  /* save the number of particles */ 
  //int nb_particles = size();
  //real t_step = time->current_step();
  //out.write(reinterpret_cast<const char*>(&nb_particles), sizeof size());
  //out.write(reinterpret_cast<const char*>(&t_step), sizeof t_step);

  boil::oout<<"save-1... " << time->current_step() << boil::endl;
  //boil::oout<<"save-1... " << time->current_step() << boil::endl -> outttt;


  std::ofstream outfile("particle_tracks.plt", std::ios_base::app);  //appendix-write-in
  //std::ofstream outfile("particle_tracks.plt");  //cover-write-in
  //outfile << " VARIABLES= xp,yp,up,vp,dp,numin" << boil::endl;

  for_p(p) {
    real dp = particles[p].d();
    real xp = particles[p].x();
    real yp = particles[p].y();
    real zp = particles[p].z(); 
    real up = particles[p].u(); 
    real vp = particles[p].v(); 
    real wp = particles[p].w();
 
    real u_old = particles[p].uvwc_old(Comp::u());  //base-velocity-?
    real v_old = particles[p].uvwc_old(Comp::v());
    real w_old = particles[p].uvwc_old(Comp::w());

    //boil::oout<<"save-d_p... " << dp << boil::endl;
    //boil::oout<<"save-xyz_p... " << xp << " " << yp << " " << zp << boil::endl;
    //boil::oout<<"save-uvw_p... " << up << " " << vp << " " << wp << boil::endl;
    //boil::oout<<"save-uvw_old... " << u_old << " " << v_old << " " << w_old << boil::endl;


    //std::ofstream outfile("particle_tracks.txt", std::ios_base::app);
    //outfile << time->current_step() << " " << xp << " " << yp << " " << up << " " << vp << " " << dp << boil::endl;
    //outfile.close();


    //outfile << " " << xp << " " << yp << " " << up << " " << vp << " " << dp << " " << (p+1) << boil::endl;
    outfile << xp << " " << yp << " " << up << " " << vp << " " << dp << " " << time->current_step() << boil::endl;

    //out.write(xp);
    //boil::oout << (particles[p].x()) -> outxp;

    /*
    out.write(reinterpret_cast<const char *> (&dp), sizeof dp);
    out.write(reinterpret_cast<const char *> (&xp), sizeof xp);
    out.write(reinterpret_cast<const char *> (&yp), sizeof yp);
    out.write(reinterpret_cast<const char *> (&zp), sizeof zp);
    out.write(reinterpret_cast<const char *> (&up), sizeof up);
    out.write(reinterpret_cast<const char *> (&vp), sizeof vp);
    out.write(reinterpret_cast<const char *> (&wp), sizeof wp);
    out.write(reinterpret_cast<const char *> (&u_old), sizeof u_old);
    out.write(reinterpret_cast<const char *> (&v_old), sizeof v_old);
    out.write(reinterpret_cast<const char *> (&w_old), sizeof w_old);*/
  }

  outfile.close();

  
  /* close a file */
  //out.close();
}

/******************************************************************************/
void Lagrangian::load(const char * nm, const int it) {

  std::string name;

  /* file name*/
  if (it == -1) 
    name = name_file(nm, ".bck", time->first_step());
  else 
    name = name_file(nm, ".bck", it);

  /* open a file */
   std::ifstream in(name.c_str(), std::ios::binary);

 /* stop if file is not present */
  if( in.rdstate() != 0 ) { 
    std::cout << "failed to open " << name << std::endl;
    std::cout << "exiting!" << std::endl;
    exit(0);
  }

  /* load the necessary data */
  int nb_particles;
  real t_step;
  in.read(reinterpret_cast<char *>(&nb_particles), sizeof size());
  in.read(reinterpret_cast<char *>(&t_step), sizeof t_step);

  *this = continuous;   
  particles.clear();

  for(int p = 0; p < nb_particles; p++) {
    real xp, yp, zp, up, vp, wp;
    real dp;
    real u_old, v_old, w_old;

    in.read(reinterpret_cast< char *> (&dp), sizeof dp);
    in.read(reinterpret_cast< char *> (&xp), sizeof xp);
    in.read(reinterpret_cast< char *> (&yp), sizeof yp);
    in.read(reinterpret_cast< char *> (&zp), sizeof zp);
    in.read(reinterpret_cast< char *> (&up), sizeof up);
    in.read(reinterpret_cast< char *> (&vp), sizeof vp);
    in.read(reinterpret_cast< char *> (&wp), sizeof wp);
    in.read(reinterpret_cast< char *> (&u_old), sizeof u_old);
    in.read(reinterpret_cast< char *> (&v_old), sizeof v_old);
    in.read(reinterpret_cast< char *> (&w_old), sizeof w_old);

    add(Particle (Position(xp,yp,zp), Diameter(dp), Position(up,vp,wp)));
  }
    
  /* close a file */
  in.close();
}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_save_load.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
