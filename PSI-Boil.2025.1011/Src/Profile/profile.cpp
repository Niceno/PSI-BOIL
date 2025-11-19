#include "profile.h"

/******************************************************************************/
Profile::Profile(const std::string name) {

  std::ifstream in;

  /*----------------+ 
  |  open the file  |
  +----------------*/
  in.open(name.c_str());

  /*----------------------------------+
  |  stop if the file is not present  |
  +----------------------------------*/
  if( in.rdstate() != 0 ) {
    boil::aout << "failed to open " << name << boil::endl;
    boil::aout << "exiting!" << boil::endl;
    exit(0);
  }

  /*--------------------+
  |  working variables  |
  +--------------------*/
  std::string item;
  std::string line;
  
  std::vector<std::string> tokens;

  /*------------------------+
  |  read the profile file  |
  +------------------------*/
  int n=0;
  while( getline(in, line) ) {

    n++;
    tokens = boil::split_string(line, " ");

    if( tokens.size() != 2 ) {
      boil::aout << "error in format of file " << name 
                 << " in line " << n << boil::endl;
      boil::aout << "exiting!" << boil::endl;
      exit(0);

    }
    
    /* transform strings into doubles */
    char * end;
    real x = strtod( tokens[0].c_str(), & end );
    real y = strtod( tokens[1].c_str(), & end );

    /* append to vectors */
    coord.push_back(x);
    value.push_back(y);
  }             

}

/******************************************************************************/
Profile::Profile(const Profile & pro) {

  coord = pro.coord;
  value = pro.value;

}

/******************************************************************************/
std::ostream & operator << (std::ostream &ost, const Profile & prof) {

  ost << boil::endl;
  for(int i=0; i<prof.coord.size(); i++) 
    ost << prof.coord[i] << " " << prof.value[i] << boil::endl;
  
  return ost;
}
