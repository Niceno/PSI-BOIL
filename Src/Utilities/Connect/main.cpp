#include "connect.h"

/* compile with: 

g++ *.cpp ../Global/libGlobal.a 

*/

/******************************************************************************/
main(int argc, char * argv[]) {

  if(argc != 4 && argc != 5) {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] 
              << " <base_name> <ext> <number_of_processor> [time_step/level]"  
              << std::endl;
    return 0;
  }

  const char * bname = argv[1];
  const int    n     = atoi(argv[3]); // number of processor

  int t = -1;
  if(argc == 5) t    = atoi(argv[4]); // time step            

  if(strcmp(argv[2], "gmv") == 0) 
    connect_gmv(bname, n, t);
  else if(strcmp(argv[2], "dat") == 0) 
    connect_tec(bname, n, t);
  else if(strcmp(argv[2], "vtk") == 0) 
    connect_vtk(bname, n, t);
  else if(strcmp(argv[2], "ib.vtk") == 0) 
    connect_ib_vtk(bname, n, t);
  else {
    std::cout << "Unsupperted format \"" << argv[2] << "\"" << std::endl;
    std::cout << "Supperted formats: \"gmv\" and \"dat\""   << std::endl;
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.6 2009/11/01 13:13:59 niceno Exp $'/
+-----------------------------------------------------------------------------*/
