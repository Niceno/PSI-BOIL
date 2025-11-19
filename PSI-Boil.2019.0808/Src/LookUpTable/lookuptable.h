/***************************************************************************//**
*  \brief A lookuptable header file for PSI-Boil                 
*
*  This class holds a table with physical properties.                  
*
*  \warning
*  
*******************************************************************************/

/* directives for the compiler */
#ifndef PROPERTYTABLE_H   
#define PROPERTYTABLE_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>                       /* first stadard C++ headers  ...     */
#include <cmath>                        /* first stadard C++ headers  ...     */
#include "../Global/global_precision.h" /* ... followed by PSI-Boil's headers */
#include "../Parallel/Out/print.h"
#include "../Ravioli/column.h"        

/////////////////////
//                 //
//  LookUpTable  //
//                 //
/////////////////////
class LookUpTable {
  public:
    //! A more elaborate constructor.
    LookUpTable(const char * filename);
    
    real look_up(const real val, const Column & col0, 
                                 const Column & col1) const;

    ~LookUpTable();

  private:
    //! Default construtor is hidden.
    LookUpTable() {}  

    std::vector < std::vector < real > > table;
};

#endif /* directive for the compiler */
