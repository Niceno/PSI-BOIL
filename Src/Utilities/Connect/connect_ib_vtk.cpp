#include "connect.h"

/******************************************************************************/
void connect_ib_vtk(const char * bname, const int n, const int t) {

  /*--------------------+
  |  working variables  |
  +--------------------*/
  std::string item;
  std::string name;
  std::string line;
  std::vector<int> nnodes_dom;
  std::vector<int> snodes_dom;
  std::vector<int> npoly_dom;
  std::vector<int> ndata_dom;
  std::vector<int> nodes;
  int nnodes=0;
  int npoly =0;
  int ndata =0;

  /*---------------+
  |  file streams  |
  +---------------*/
  std::ofstream out;
  std::ifstream * in = new std::ifstream[n];
  nnodes_dom.resize(n);
  snodes_dom.resize(n);
  npoly_dom.resize(n);
  ndata_dom.resize(n);
  nnodes=0;

  /*---------------------------------------+
  |  >> browse through files to open them  |
  +---------------------------------------*/
  for(int f=0; f<n; f++) {
    std::string current = name_file(bname, ".ib.vtk", t, f);
    std::cout << "reading: " << current << std::endl;

    in[f].open(current.c_str());

    /* stop if file is not present */
    if( in[f].rdstate() != 0 ) {
      std::cout << "failed to open " << current << std::endl;
      std::cout << "exiting!" << std::endl;
      exit(0);
    }
  }

  /*----------------------+
  |  << start the output  |
  +----------------------*/
  std::string out_name = name_file(bname, ".ib.vtk", t);
  std::cout << "creating: " << out_name << std::endl;
  out.open( out_name.c_str() );

  out << "# vtk DataFile Version 2.0"      << std::endl;
  out << "this is a comment use it wisely" << std::endl;
  out << "ASCII"                           << std::endl;
  out << "DATASET POLYDATA"                << std::endl;

  /*--------------------------------+
  |  >> read total number of nodes  | 
  +--------------------------------*/
  for(int f=0; f<n; f++) {
    getline(in[f], line);         /* "# vtk DataFile Version 2.0" */
    getline(in[f], line);         /* comment                      */
    getline(in[f], line);         /* "ASCII"                      */
    getline(in[f], line);         /* "DATASET POLYDATA"           */

    in[f] >> item;                /* "POINTS"        */
    in[f] >> nnodes_dom[f];       /* number of nodes */
    nnodes += nnodes_dom[f];      /* total number of nodes */
    getline(in[f], line);         /* finish the line */
  }
  out << "POINTS " << nnodes << " float" << std::endl;
  
  /*-------------------------------------+ 
  |  get starting nodes for each domain  |
  +-------------------------------------*/
  snodes_dom[0]=0;
  for(int f=1; f<n; f++) {
    snodes_dom[f]=snodes_dom[f-1]+nnodes_dom[f-1];
  }

  /*--------------------+
  |  << >> coordinates  |
  +--------------------*/
  for(int f=0; f<n; f++) {
    for(int i=0; i<nnodes_dom[f]; i++) {
      getline(in[f], line);   
      out << line << std::endl; 
    }    
  }

  /*-----------------------------------+
  |  >> read total number of polygons  |
  +-----------------------------------*/
  npoly = 0;
  ndata = 0;
  for(int f=0; f<n; f++) {
    in[f] >> item;                /* "POLYGONS"      */
    in[f] >> npoly_dom[f];        /* number of polygons */
    in[f] >> ndata_dom[f];        /* number of datums   */
    npoly += npoly_dom[f];        /* number of polygons */
    ndata += ndata_dom[f];        /* number of datums   */
    getline(in[f], line);         /* finish the line */
  }
  out << "POLYGONS " << npoly << " " << ndata << std::endl;

  /*-----------------+
  |  << >> polygons  |
  +-----------------*/
  nodes.resize(128);
  for(int f=0; f<n; f++) {

    for(int p=0; p<npoly_dom[f]; p++) {

      /* >> */
      int nn; in[f] >> nn;
      for(int i=0; i<nn; i++) in[f] >> nodes[i];
      getline(in[f], line); /* finish the line */
 
      /* << */
      out << nn << " ";          
      for(int i=0; i<nn; i++) out << nodes[i] + snodes_dom[f] << " ";
      out << std::endl;
    }
  }
}
