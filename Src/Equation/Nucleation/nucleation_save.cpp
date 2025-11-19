#include "nucleation.h"

/******************************************************************************/
void Nucleation::save(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  std::ofstream out(name.c_str(), std::ios::binary);
  
  /* save the necessary data */
  save(out);
  
  /* close a file */
  out.close();
}

/******************************************************************************/
void Nucleation::save(std::ofstream & out) {

  int sites_size = size();
  out.write(reinterpret_cast<const char *> (&sites_size), sizeof(int));

  for(int ns=0; ns<size(); ns++){
    real x = sites[ns].x();
    real y = sites[ns].y();
    real z = sites[ns].z();
    real c1 = sites[ns].active_tpr();
    real c2 = sites[ns].zplant();
    real tseed = sites[ns].time_plant_clr();
    bool bseedp = sites[ns].seed_prev();
    bool act    = sites[ns].active();
    out.write(reinterpret_cast<const char *> (&x), sizeof(real));
    out.write(reinterpret_cast<const char *> (&y), sizeof(real));
    out.write(reinterpret_cast<const char *> (&z), sizeof(real));
    out.write(reinterpret_cast<const char *> (&c1), sizeof(real));
    out.write(reinterpret_cast<const char *> (&c2), sizeof(real));
    out.write(reinterpret_cast<const char *> (&tseed), sizeof(real));
    out.write(reinterpret_cast<const char *> (&bseedp), sizeof(bool));
    out.write(reinterpret_cast<const char *> (&act), sizeof(bool));
    if (pre_heat_sink()) {
      real sqe = sites[ns].sum_sink_energy();
      real tTact = sites[ns].time_Tact();
      bool bqsink = sites[ns].qsink();
      out.write(reinterpret_cast<const char *> (&sqe),    sizeof(real));
      out.write(reinterpret_cast<const char *> (&tTact),  sizeof(real));
      out.write(reinterpret_cast<const char *> (&bqsink), sizeof(bool));
    }
  }

}
