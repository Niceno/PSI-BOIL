#include "body.h"
#include "../Global/global_approx.h"

/******************************************************************************/
Body::Body(const std::string nm) {

  name = nm;

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
  
  real xp[3], yp[3], zp[3], np[3];

  getline(in, line);                         /* header          */

  for(;;) {

    in >> item;                              /* "facet" or "endsolid" */
    if(item == "endsolid") break;  
    in >> item;                              /* "normal"  */

      in >> np[0]; in >> np[1]; in >> np[2]; /* n_x, n_y, n_z   */
      getline(in, line);                     /* finish the line */

      getline(in, line);                     /* "outer loop"    */
        in >> item; in >> xp[0]; in >> yp[0]; in >> zp[0]; getline(in, line); 
        in >> item; in >> xp[1]; in >> yp[1]; in >> zp[1]; getline(in, line); 
        in >> item; in >> xp[2]; in >> yp[2]; in >> zp[2]; getline(in, line); 
      getline(in, line);                     /* "endloop"       */

      Polygon * t = new Polygon(3, xp, yp, zp/*, np*/); // CHECKING
      polys.push_back(*t);
      delete t;

    getline(in, line);                       /* "endfacet"      */

  }

  /*----------------------+
  |  initialize pointers  |
  +----------------------*/
  sca = NULL;
  vec = NULL;
  bdist  = NULL;
  stmp   = NULL;
  dflag  = NULL;
  ux = NULL;
  uy = NULL;
  uz = NULL;
  vecoff = NULL;

  //debug: boil::aout << (int)polys.size() << boil::endl;
}
