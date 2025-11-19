#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::cell_init() {
/*---------------------------------------------------------------------+
|  this create fictitious coarse cells in the domain.                  |
|  the only reason that we implement it is to speedup                  |
|  the collision algorithm.                                            |
|  for a given particle, we only check if it will collide              |  
|  with the particles surrounding it, instead of checking              |
|  with all the particles in the domain.                               |
|  for more details, check Hassan Badreddine's thesis section (2.3.7)  | 
+---------------------------------------------------------------------*/

  /* setup the fictitious domain */
  const real gx = dom->global_max_x() - dom->global_min_x();
  const real gy = dom->global_max_y() - dom->global_min_y();
  const real gz = dom->global_max_z() - dom->global_min_z();

  NX_coarse = int (gx/list_diameter);
  NY_coarse = int (gy/list_diameter);
  NZ_coarse = int (gz/list_diameter);

  diam_x  = gx / double (NX_coarse);
  diam_y  = gy / double (NY_coarse);
  diam_z  = gz / double (NZ_coarse);
  OPR(diam_x); OPR(diam_y); OPR(diam_z); 

  if(NX_coarse == 0) {
    boil::oout << "NX_coarse is 0..change list_diameter in lagrangian.cpp"
               << '\n';
    exit(0);
  }
  if(NY_coarse == 0) { 
    boil::oout << "NY_coarse is 0.. change list_diameter in lagrangian.cpp" 
               << '\n';
    exit(0);
  }
  if(NZ_coarse == 0) { 
    boil::oout << "NZ_coarse is 0.. change list_diameter in lagrangian.cpp" 
               << '\n';
    exit(0);
  }

  /* declare the cell matrix */
  cell = new int ** [NX_coarse];
  for (int i = 0; i < NX_coarse; i++) {
    cell[i] = new int * [NY_coarse];
    for(int j = 0; j <  NY_coarse; j++) cell[i][j] = new int [NZ_coarse];
  }

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_cell_init.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
