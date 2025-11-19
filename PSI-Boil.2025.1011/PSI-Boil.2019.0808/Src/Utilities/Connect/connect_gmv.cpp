#include "connect.h"

/******************************************************************************/
void connect_gmv(const char * bname, const int n, const int t) {

  /*--------------------+
  |  working variables  |
  +--------------------*/
  std::string item, item_1, item_2;
  std::string line;
  std::vector<int> nnodes_dom;       /* number of nodes in a domain */
  std::vector<int> ncells_dom;       /* number of cells in a domain */
  std::vector<int> first_node_dom;   /* number of cells in a domain */
  std::vector<int> cell_node;        /* cell nodes */
  nnodes_dom.resize(n);
  ncells_dom.resize(n);
  first_node_dom.resize(n);
  cell_node.resize(8);
  int nnodes = 0;
  int ncells = 0;
  std::vector<Node> node_array; /* array of nodes */
  std::vector<int> new_node_number;

  /*---------------+
  |  file streams  |
  +---------------*/
  std::ofstream out;
  std::ifstream * in = new std::ifstream[n];

  /*--------------------------+
  |  browse through files to  |
  |   read node coordinates   |
  |     and cell numbers      |
  +--------------------------*/
  nnodes = 0; /* count nodes */
  for(int f=0; f<n; f++) {
    std::string current = name_file(bname, ".gmv", t, f);
    std::cout << "reading: " << current << std::endl;

    in[f].open(current.c_str());

    /* stop if file is not present */
    if( in[f].rdstate() != 0 ) {
      std::cout << "failed to open " << current << std::endl;
      std::cout << "exiting!" << std::endl;
      exit(0);
    }

    /*----------------------+
    |  >> node coordinates  | 
    +----------------------*/
    getline(in[f], line); /* "gmvimput ascii" */
    in[f] >> item;        /* "nodev"  */
    getline(in[f], item); /* number */

    nnodes_dom[f] = atoi(item.c_str());
    std::cout << nnodes_dom[f] << " nodes" << std::endl;

    real x,y,z;

    for(int l=0; l<nnodes_dom[f]; l++) { /* read all coordinates */
      nnodes++;                    /* one more node */
      in[f] >> x >> y >> z; 
    
      node_array.push_back( Node(nnodes, x, y, z) );
    }

    /*------------------+
    |  >> cell numbers  | 
    +------------------*/
    in[f] >> item;        /* "cells"  */
    getline(in[f], item); /* number */

    ncells_dom[f] = atoi(item.c_str()); /* get cell number */
    ncells += ncells_dom[f];
    std::cout << ncells_dom[f] << " cells" << std::endl;
  }

  /*---------------------+
  |                      |
  |  renumber the nodes  |
  |                      |
  +---------------------*/
  new_node_number.resize(nnodes+1);
  int nnodes_old = nnodes;

  stable_sort(node_array.begin(), node_array.end());

  nnodes=1;
  node_array[0].new_number(nnodes);
  for(int i=1; i<node_array.size(); i++) {
    if(node_array[i-1] < node_array[i]) nnodes++;
    node_array[i].new_number(nnodes);
  }
  std::cout << nnodes << std::endl;

  for(int i=0; i<node_array.size(); i++) {
    int oldn = node_array[i].old_number();
    int newn = node_array[i].new_number();
    new_node_number[oldn] = newn;
  }

  /*--------------------------------------+
  |  estimate first node for each domain  |
  +--------------------------------------*/
  first_node_dom[0] = 0;
  for(int f=1; f<n; f++)
    first_node_dom[f] = first_node_dom[f-1] + nnodes_dom[f-1];

  /*-------------------+
  |  start the output  |
  +-------------------*/
  std::string out_name = name_file(bname, ".gmv", t);
  std::cout << "creating: " << out_name << std::endl;
  out.open( out_name.c_str() );

  /*-----------+
  |  << nodes  |
  +-----------*/
  out << "gmvinput ascii"   << std::endl; 
  out << "nodev " << nnodes << std::endl; 

  /* <- write all coordinates */
  out << node_array[0] << std::endl;
  for(int i=1; i<node_array.size(); i++) 
    if(node_array[i-1] < node_array[i]) 
      out << node_array[i] << std::endl;

  /*--------------+
  |  >> << cells  |
  +--------------*/
  out << "cells " << ncells << std::endl;
  for(int f=0; f<n; f++) {

    for(int l=0; l<ncells_dom[f]; l++) { /* write all coordinates */
      in[f] >> item_1;
      in[f] >> item_2;
      getline(in[f], line);      
      if(item_1 == "hex") { 
        out << item_1 << " " << item_2 << std::endl;  /* <- hex 8 */
        for(int k=0; k<8; k++) in[f] >> cell_node[k];  
        for(int k=0; k<8; k++) 
          out << new_node_number[ cell_node[k]+first_node_dom[f] ] << " "; 
        out << std::endl;           
        getline(in[f], line); /* finish the line */
      }
      else if(item_1 == "general") { 
        out << item_1 << " " << item_2 << std::endl;  /* <- general 1 */
        int n;
        in[f] >> n;                            
        getline(in[f], line); /* finish the line */
        out << n << std::endl;  /* <- n  */
        for(int k=0; k<n; k++) in[f] >> cell_node[k];  
        for(int k=0; k<n; k++) 
          out << new_node_number[ cell_node[k]+first_node_dom[f] ] << " "; 
        out << std::endl;           
        getline(in[f], line); /* finish the line */
      }
    }
  }

  /*-----------------------------------------+
  |  >> << materials, variables, velocities  |
  +-----------------------------------------*/
  for(;;) {

    bool mate = false;
    bool vars = false;
    bool velo = false;
    bool end  = false;
    for(int i=0; i<n; i++) {
      in[i] >> item;        /* -> variables / velocities  */
      if(item == "material")  mate = true;
      if(item == "variables") vars = true;
      if(item == "velocity")  velo = true;
      if(item == "endgmv")    end  = true;
    }

    /*------+
    |  end  |
    +------*/
    if(end) {
      for(int i=0; i<n; i++) in[i].close();
      out << "endgmv" << std::endl;
      out.close();
      std::cout << "done!" << std::endl;
      return; 
    }

    std::cout << "processing: " << item << std::endl;
  
    /*------------+
    |  materials  |
    +------------*/
    if(mate) {
      int nm, pos;
      for(int f=0; f<n; f++) {
        in[f] >> nm;          /* -> number of materials */
        in[f] >> pos;         /* -> position (0) */
        getline(in[f], line); /* finish the line */
  
        if(f==0) out << item << " "        /* <- material */
                     << nm << " "          /* <- number of materials */
                     << pos << std::endl;  /* <- position */
  
        for(int m=0; m<nm; m++) {
          getline(in[f], line);                /* -> material name */
          if(f==0) {out << line << std::endl;  /* <- material name */
                    std::cout << line << std::endl;} 
  
        }
  
        for(int l=0; l<ncells_dom[f]; l++) { /* write all material tags */
          getline(in[f], line);      /* -> value */
          out << line << std::endl;  /* <- value */
        }
      }
    }

    /*------------+
    |  variables  |
    +------------*/
    if(vars) {
      for(int f=0; f<n; f++) {
        getline(in[f], line); /* finish the line */
        if(f==0) out << item << std::endl;  /* <- "variables" */
      }
  
      for(;;) {
        for(int f=0; f<n; f++) {
          in[f] >> item_1;                       /* -> name (or endvars) */
          if(item_1 == "endvars") goto endloop;
          in[f] >> item_2;                       /* -> variable position */
          getline(in[f], line); /* finish the line */
          if(f==0) 
            out << item_1 << " " << item_2 << std::endl; /* <- name, pos */

          for(int l=0; l<ncells_dom[f]; l++) { /* write all variables */
            getline(in[f], line);      /* -> value */
            out << line << std::endl;  /* <- value */
          }
        }
      }
endloop:
      out << "endvars" << std::endl;
    }

    /*-----------+
    |  velocity  |
    +-----------*/
    if(velo) {
      for(int f=0; f<n; f++) {
        getline(in[f], line); /* finish the line */
        if(f==0) out << item << line << std::endl;  /* <- "velocity 0" */
      }

      for(int m=0; m<3; m++) { /* through velocity components */
        for(int f=0; f<n; f++) {
          for(int l=0; l<ncells_dom[f]; l++) { /* write all components */
            getline(in[f], line);      /* -> value */
            out << line << std::endl;  /* <- value */
          }
        }
      }
    }
  }

}
