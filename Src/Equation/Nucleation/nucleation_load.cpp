#include "nucleation.h"
using namespace std;

/******************************************************************************/
void Nucleation::load(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  std::ifstream in(name.c_str(), std::ios::binary);
  
  /* stop if file is not present */
  if( in.rdstate() != 0 ) {
    std::cout << "failed to open " << name << std::endl;
    std::cout << "exiting!" << std::endl;
    exit(0);
  }

  /* load the necessary data */
  load(in);
  
  /* close a file */
  in.close();
}

/******************************************************************************/
void Nucleation::load(std::ifstream & in) {

  int site_size;
  real x_saved, y_saved, z_saved, c1_saved, c2_saved, t_saved;
  bool bseedp, act;
  
  in.read(reinterpret_cast<char *> (&site_size), sizeof(int));
  boil::oout<<"nucleation_load: total number of seed = "<<site_size<<"\n";

  for(int ns=0; ns<site_size; ns++){
    in.read(reinterpret_cast<char *> (&x_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&y_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&z_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&c1_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&c2_saved), sizeof(real));
    add(Site(x_saved, y_saved, z_saved, c1_saved, c2_saved));
    in.read(reinterpret_cast<char *> (&t_saved), sizeof(real));
    sites[ns].set_time_plant_clr(t_saved);
    in.read(reinterpret_cast<char *> (&bseedp), sizeof(bool));
    sites[ns].set_seed_prev(bseedp);
    in.read(reinterpret_cast<char *> (&act), sizeof(bool));
    sites[ns].set_active(act);
    if (pre_heat_sink()) {
      real sqe, tTact;
      bool bqsink;
      in.read(reinterpret_cast<char *> (&sqe),    sizeof(real));
      in.read(reinterpret_cast<char *> (&tTact),  sizeof(real));
      in.read(reinterpret_cast<char *> (&bqsink), sizeof(bool));
      sites[ns].set_sum_sink_energy(sqe);
      sites[ns].set_time_Tact(tTact);
      sites[ns].set_qsink(bqsink);
    }

    if (ns%1000==1) {
      boil::oout<<"nucleation_load: "<<ns<<" "<<x_saved<<" "<<y_saved<<" "
                <<z_saved<<" "<<c1_saved<<" "<<c2_saved<<" "
                <<t_saved<<" "<<bseedp<<" "<<act<<"\n";
    }
  }

}
