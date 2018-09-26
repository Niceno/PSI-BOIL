#include "lookuptable.h" /* class header must be included */

/***************************************************************************//**
*  This is a non-trivial constructor. Note that this comment is written in
*  Doxygen style and therefore in a grammatically correct English language.            
*******************************************************************************/
LookUpTable :: LookUpTable(const char * filename) {

  std::ifstream in;
  in.open(filename);

  /* stop if file is not present */
  if( in.rdstate() != 0 ) {
    boil::oout << "failed to open: " << filename << ". exiting!" << boil::endl;
    exit(0);
  }

  std::string line;

  /* skip first two lines */
  getline(in, line);
  getline(in, line);

  /* read the rest of the lines */
  while ( getline(in, line) ) {
    std::vector < real > data;    /* row of data                             */
    real value;                   /* datum                                   */
    std::istringstream iss(line); /* manipulate string as if it was stream   */
    while (iss >> value)          /* read all entries in the current row ... */ 
      data.push_back(value);      /* ... and push them into "data"           */
    table.push_back(data);        /* add the entire row into table           */
  }

  // boil::oout << "table size: " << table.size() << boil::endl;
}

/******************************************************************************/
LookUpTable :: ~LookUpTable() {
/*--------------------------------------------------------------------------+
|  destructor is also placed here.                                          |
|  this is a psi-boil-style comment, written with lower-case letters only.  |
+--------------------------------------------------------------------------*/

}

/*-----------------------------------------------------------------------------+
 '$Id: lookuptable.cpp,v 1.1 2011/05/28 14:30:01 niceno Exp $'/ 
+-----------------------------------------------------------------------------*/
