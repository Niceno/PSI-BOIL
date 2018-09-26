#include "floodfill.h"

/******************************************************************************/
void Floodfill::save(const char * nm, const int it) {

  std::string name = name_file(nm, ".bck", it);

  std::ofstream out(name.c_str(), std::ios::binary);

  int nb_regions = m_vectrgns.size();
  real t_step = time->current_step();
  int sz_idavail = v_idavail.size();
  out.write(reinterpret_cast<const char *>(&nb_regions), sizeof nb_regions);
  out.write(reinterpret_cast<const char *>(&t_step), sizeof t_step);
  out.write(reinterpret_cast<const char *>(&prgn_inc), sizeof prgn_inc);
  out.write(reinterpret_cast<const char *>(&nrgn_inc), sizeof nrgn_inc);
  out.write(reinterpret_cast<const char *>(&sz_idavail), sizeof sz_idavail);
  out.write(reinterpret_cast<const char *>(&oldtime), sizeof oldtime);

  for (int i=0; i<v_idavail.size(); i++) {
    int a = v_idavail[i];
    out.write(reinterpret_cast<const char *> (&a), sizeof a);
  }
    
  for (int i=0; i<m_vectrgns.size(); i++) { //browse through all regions in m_vectrgns
    real rx = m_vectrgns[i].x();
    real ry = m_vectrgns[i].y();
    real rz = m_vectrgns[i].z();
    real ox = m_vectrgns[i].xold();
    real oy = m_vectrgns[i].yold();
    real oz = m_vectrgns[i].zold();
    real ru = m_vectrgns[i].u();
    real rv = m_vectrgns[i].v();
    real rw = m_vectrgns[i].w();
    int rcellvol = m_vectrgns[i].cellvol();
    int rid = m_vectrgns[i].id();
    int rhiding = m_vectrgns[i].hiding();
    int rtsteps_hidden = m_vectrgns[i].get_tsteps_hidden();

    out.write(reinterpret_cast<const char *> (&rx), sizeof rx);
    out.write(reinterpret_cast<const char *> (&ry), sizeof ry);
    out.write(reinterpret_cast<const char *> (&rz), sizeof rz);
    out.write(reinterpret_cast<const char *> (&ox), sizeof ox);
    out.write(reinterpret_cast<const char *> (&oy), sizeof oy);
    out.write(reinterpret_cast<const char *> (&oz), sizeof oz);
    out.write(reinterpret_cast<const char *> (&ru), sizeof ru);
    out.write(reinterpret_cast<const char *> (&rv), sizeof rv);
    out.write(reinterpret_cast<const char *> (&rw), sizeof rw);
    out.write(reinterpret_cast<const char *> (&rcellvol), sizeof rcellvol);
    out.write(reinterpret_cast<const char *> (&rid), sizeof rid);
    out.write(reinterpret_cast<const char *> (&rhiding), sizeof rhiding);
    out.write(reinterpret_cast<const char *> (&rtsteps_hidden), sizeof rtsteps_hidden);
  }
  
  /* close a file */
  out.close();
}

/******************************************************************************/
void Floodfill::load(const char * nm, const int it) {

  std::string name;

  if (it <= -1)  
    name = name_file(nm, ".bck", time->first_step());
  else 
    name = name_file(nm, ".bck", it);

   std::ifstream in(name.c_str(), std::ios::binary);

  if( in.rdstate() != 0 ) { 
    std::cout << "failed to open " << name << std::endl;
    std::cout << "exiting!" << std::endl;
    exit(0);
  }

  int nb_regions;
  real t_step;
  int sz_idavail;
  in.read(reinterpret_cast<char *>(&nb_regions), sizeof nb_regions);
  in.read(reinterpret_cast<char *>(&t_step), sizeof t_step);
  in.read(reinterpret_cast<char *>(&prgn_inc), sizeof prgn_inc);
  in.read(reinterpret_cast<char *>(&nrgn_inc), sizeof nrgn_inc); 
  in.read(reinterpret_cast<char *>(&sz_idavail), sizeof sz_idavail);
  in.read(reinterpret_cast<char *>(&oldtime), sizeof oldtime);

  v_idavail.clear();
  for (int i=0; i<sz_idavail; i++) {
    int a;
    in.read(reinterpret_cast< char *> (&a), sizeof a);
    v_idavail.push_back(a);
  }
    
  m_vectrgns.clear(); 
  for (int i=0; i<nb_regions; i++) { 
    real rx, ry, rz, ox, oy, oz, ru, rv, rw;
    int rcellvol, rid, rtsteps_hidden;
    int rhiding;

    in.read(reinterpret_cast< char *> (&rx), sizeof rx);
    in.read(reinterpret_cast< char *> (&ry), sizeof ry);
    in.read(reinterpret_cast< char *> (&rz), sizeof rz);
    in.read(reinterpret_cast< char *> (&ox), sizeof ox);
    in.read(reinterpret_cast< char *> (&oy), sizeof oy);
    in.read(reinterpret_cast< char *> (&oz), sizeof oz);
    in.read(reinterpret_cast< char *> (&ru), sizeof ru);
    in.read(reinterpret_cast< char *> (&rv), sizeof rv);
    in.read(reinterpret_cast< char *> (&rw), sizeof rw);
    in.read(reinterpret_cast< char *> (&rcellvol), sizeof rcellvol);
    in.read(reinterpret_cast< char *> (&rid), sizeof rid);
    in.read(reinterpret_cast< char *> (&rtsteps_hidden), sizeof rtsteps_hidden);
    in.read(reinterpret_cast< char *> (&rhiding), sizeof rhiding);

    Region & added_rgn = getregion(rid); //sets id
    added_rgn.xyz(rx, ry, rz);  
    added_rgn.oldxyz(ox, oy, oz);  
    added_rgn.uvw(ru, rv, rw);
    added_rgn.cellvol(rcellvol);
    added_rgn.set_nts_hidden(rtsteps_hidden);
    added_rgn.hiding(rhiding);

  }
    
  in.close();
}

/*-----------------------------------------------------------------------------+
 '$Id: floodfill_save_load.cpp,v 1.1 2018/02/16 19:07:11 sato Exp $'/
+-----------------------------------------------------------------------------*/
