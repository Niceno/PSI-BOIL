#include "custom.h"

namespace boil {
  std::string save_backup(const int ts, const bool irregular,
                   const Times & time,
                   const std::vector<Scalar*> & scalars,
                   const std::vector<std::string> & scalar_names,
                   const std::vector<Vector*> & vectors,
                   const std::vector<std::string> & vector_names,
                   const std::vector<Nucleation*> & nucls,
                   const std::vector<std::string> & nucl_names,
                   const std::vector<CIPCSL2*> & cipcsl2s,
                   const std::vector<std::string> & cipcsl2_names,
                   const std::vector<real*> & store_values) {
    /* file name */
    std::stringstream ss;
    if(!irregular) {
      ss<<"time-"<<ts<<".txt";
    } else {
      ss<<"time.txt";
    }   
    std::string fname = ss.str();

    /* save scalars */
    for(int i(0); i<scalars.size(); ++i) {
       scalars[i]->save(scalar_names[i].c_str(),ts);
    }

    /* save vectors */
    for(int i(0); i<vectors.size(); ++i) {
       vectors[i]->save(vector_names[i].c_str(),ts);
    }

    /* save nucleations */
    for(int i(0); i<nucls.size(); ++i) {
       nucls[i]->save(nucl_names[i].c_str(),ts);
    }

    /* save CIPCSL2s */
    for(int i(0); i<cipcsl2s.size(); ++i) {
       cipcsl2s[i]->save(cipcsl2_names[i].c_str(),ts);
    }

    /* generate time file */
    if( boil::cart.iam()==0) {
      std::fstream output;
      output << std::setprecision(16);
      output.open(fname, std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;

      /* save tracked values */
      for(auto val : store_values)
        output << *val << boil::endl;
      output.close();
    }

    return fname;
  }

  bool load_backup(const std::string & fname,
                   int & ts, Times & time,
                   const std::vector<Scalar*> & scalars,
                   const std::vector<std::string> & scalar_names,
                   const std::vector<Vector*> & vectors,
                   const std::vector<std::string> & vector_names,
                   const std::vector<Nucleation*> & nucls,
                   const std::vector<std::string> & nucl_names,
                   const std::vector<CIPCSL2*> & cipcsl2s,
                   const std::vector<std::string> & cipcsl2_names,
                   const std::vector<real*> & store_values) {
    ts = 0;
    std::fstream input;
    input.open(fname,std::ios::in);
    bool open_failed = input.fail();

    if(!open_failed) {
      /* simulation time */
      real t,dtf;
      input >> ts;
      input >> t;
      input >> dtf;
      time.first_step(ts);
      time.current_time(t);
      time.set_dt(dtf);

      /* load tracked values */
      //if(store_values != NULL)
        for(auto val : store_values)
          input >> *val;

      /* load scalars */
      for(int i(0); i<scalars.size(); ++i) {
         scalars[i]->load(scalar_names[i].c_str(),ts);
      }

      /* load vectors */
      for(int i(0); i<vectors.size(); ++i) {
         vectors[i]->load(vector_names[i].c_str(),ts);
      }
    
      /* load nucleations */
      for(int i(0); i<nucls.size(); ++i) {
         nucls[i]->load(nucl_names[i].c_str(),ts);
      }

      /* load CIPCSL2s */
      for(int i(0); i<cipcsl2s.size(); ++i) {
         cipcsl2s[i]->load(cipcsl2_names[i].c_str(),ts);
      }

    }
    input.close();

    return (!open_failed);
  }

  void rm_backup(const int ts, 
                 const std::vector<Scalar*> & scalars,
                 const std::vector<std::string> & scalar_names,
                 const std::vector<Vector*> & vectors,
                 const std::vector<std::string> & vector_names,
                 const std::vector<Nucleation*> & nucls,
                 const std::vector<std::string> & nucl_names,
                 const std::vector<CIPCSL2*> & cipcsl2s,
                 const std::vector<std::string> & cipcsl2_names,
                 const std::vector<real*> & values) {
    /* rm scalars */
    for(int i(0); i<scalars.size(); ++i) {
       scalars[i]->rm(scalar_names[i].c_str(),ts);
    }

    /* rm vectors */
    for(int i(0); i<vectors.size(); ++i) {
       vectors[i]->rm(vector_names[i].c_str(),ts);
    }

    /* rm nucleations */
    for(int i(0); i<nucls.size(); ++i) {
       nucls[i]->rm(nucl_names[i].c_str(),ts);
    }

    /* rm CIPCSL2s */
    for(int i(0); i<cipcsl2s.size(); ++i) {
       cipcsl2s[i]->rm(cipcsl2_names[i].c_str(),ts);
    }

    return;
  }

} /* namespace */
