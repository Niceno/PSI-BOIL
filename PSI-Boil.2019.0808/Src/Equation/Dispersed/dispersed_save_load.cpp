#include "dispersed.h"

/******************************************************************************/
void Dispersed::save(const char * nm, const int it) {

  /* file name*/
  std::string name = name_file(nm, ".bck", it);

  /* open a file */
  std::ofstream out(name.c_str(), std::ios::binary);

  /* save the number of particles */ 
  int nb_particles = size();
  real t_step = time->current_step();
  out.write(reinterpret_cast<const char*>(&nb_particles), sizeof size());
  out.write(reinterpret_cast<const char*>(&t_step), sizeof t_step);

  for_p(p) {
    real dp = particles[p].d();
    real xp = particles[p].x();
    real yp = particles[p].y();
    real zp = particles[p].z(); 
    real up = particles[p].u(); 
    real vp = particles[p].v(); 
    real wp = particles[p].w();
 
    real u_old = particles[p].uvwc_old(Comp::u());
    real v_old = particles[p].uvwc_old(Comp::v());
    real w_old = particles[p].uvwc_old(Comp::w());

    out.write(reinterpret_cast<const char *> (&dp), sizeof dp);
    out.write(reinterpret_cast<const char *> (&xp), sizeof xp);
    out.write(reinterpret_cast<const char *> (&yp), sizeof yp);
    out.write(reinterpret_cast<const char *> (&zp), sizeof zp);
    out.write(reinterpret_cast<const char *> (&up), sizeof up);
    out.write(reinterpret_cast<const char *> (&vp), sizeof vp);
    out.write(reinterpret_cast<const char *> (&wp), sizeof wp);
    out.write(reinterpret_cast<const char *> (&u_old), sizeof u_old);
    out.write(reinterpret_cast<const char *> (&v_old), sizeof v_old);
    out.write(reinterpret_cast<const char *> (&w_old), sizeof w_old);
  }
  
  /* close a file */
  out.close();
}

/******************************************************************************/
void Dispersed::load(const char * nm, const int it) {

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
