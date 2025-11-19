#include "connect.h"

/******************************************************************************/
void connect_tec(const char * bname, const int n, const int t) {

  /*--------------------+
  |  working variables  |
  +--------------------*/
  std::string item;
  std::string type;
  std::string line;
  std::string line_var;
  std::string line_com;
  std::string store_1, store_2;
  std::string format;
  std::vector<real> xyz;

  /*---------------+
  |  file streams  |
  +---------------*/
  std::ofstream out;
  std::ifstream * in = new std::ifstream[n];

  /*---------------------------------------+
  |  >> browse through files to open them  |
  +---------------------------------------*/
  for(int f=0; f<n; f++) {
    std::string current = name_file(bname, ".dat", t, f);
    std::cout << "reading: " << current << std::endl;

    in[f].open(current.c_str());

    /* stop if file is not present */
    if( in[f].rdstate() != 0 ) {
      std::cout << "failed to open " << current << std::endl;
      std::cout << "exiting!" << std::endl;
      exit(0);
    }
  }

  /*---------------------------------------------+
  |  >> browse through files to get i,j,k range  | 
  +---------------------------------------------*/
  std::vector<int>  is, ie, js, je, ks, ke;
  is.resize(n); ie.resize(n);
  js.resize(n); je.resize(n);
  ks.resize(n); ke.resize(n);
  int gis=INT_MAX, gie=-INT_MAX, 
      gjs=INT_MAX, gje=-INT_MAX, 
      gks=INT_MAX, gke=-INT_MAX;
  int nodal=0;

  for(int f=0; f<n; f++) {
    getline(in[f], line);         /* "# FILE CREATED WITH PSI-BOIL" */

    /* get format */
    in[f] >> item; in[f] >> format; getline(in[f], line);
    if( format == "NODAL" ) nodal=1;

    /* get ranges */
    in[f] >> item; in[f] >> item; in[f] >> is[f]; in[f] >> ie[f];
    in[f] >> item; in[f] >> item; in[f] >> js[f]; in[f] >> je[f];
    in[f] >> item; in[f] >> item; in[f] >> ks[f]; in[f] >> ke[f];

    if(ie[f]>=is[f] && je[f]>=js[f] && ke[f]>=ks[f]) {
      if(is[f] < gis) gis=is[f]; if(ie[f] > gie) gie=ie[f];
      if(js[f] < gjs) gjs=js[f]; if(je[f] > gje) gje=je[f];
      if(ks[f] < gks) gks=ks[f]; if(ke[f] > gke) gke=ke[f];
    }

    getline(in[f], line);         /* finish the line */
  }

  /*-------------------+
  |  >> variable info  | 
  +-------------------*/
  int n_var;
  for(int f=0; f<n; f++) {
    getline(in[f], line_var);     /* " VARIABLES=... */
 
    /* count variables (by counting the quotations) */
    n_var = 0;
    for(int i=0; i < line_var.length(); i++)
      if( line_var.at(i) == '\"' ) n_var++;
    n_var /= 2;
    n_var -= 3;
    std::cout << "number of variables (excluding coordinates): " 
              << n_var << std::endl;  

    in[f] >> item;                /* " ZONE" */
    in[f] >> item; in[f] >> item; in[f] >> item; 
    in[f] >> store_1; 
    if(n_var > 0)
      in[f] >> store_2;
    getline(in[f], line);         /* " finish the line */
  }

  /*----------------------+
  |  << start the output  |
  +----------------------*/
  std::string out_name = name_file(bname, ".dat", t);
  std::cout << "creating: " << out_name << std::endl;
  out.open( out_name.c_str() );

  out << "# FILE CREATED WITH PSI-BOIL"     << std::endl;
  out << "# These three lines are  " << std::endl;
  out << "# for compatibility with " << std::endl;
  out << "# MovieScript utility.   " << std::endl;
  if(nodal) out << "# NODAL FORMAT "        << std::endl;
  else      out << "# CELLCENTERED FORMAT " << std::endl;

  out << line_var << std::endl;

  out << "  ZONE I=" << gie+1-nodal   -gis+1
      <<      ", J=" << gje+1-nodal   -gjs+1
      <<      ", K=" << gke+1-nodal   -gks+1
      << " " << store_1 << " " << store_2 << std::endl;

  /*-------------------------+
  |  >> << node coordinates  | 
  +-------------------------*/
  for(int m=0; m<3; m++) {
    real coord;

    if(m==0) xyz.resize(gie+1);
    if(m==1) xyz.resize(gje+1);
    if(m==2) xyz.resize(gke+1);

    /* >> */
    for(int f=0; f<n; f++) {
      if(ie[f]>=is[f] && je[f]>=js[f] && ke[f]>=ks[f]) {
        getline(in[f], line_com);     /* "# *-COORDINATE" */
        for(int k=ks[f]-1; k<=ke[f]-nodal; k++) 
          for(int j=js[f]-1; j<=je[f]-nodal; j++)
            for(int i=is[f]-1; i<=ie[f]-nodal; i++) {
              if(m==0) in[f] >> xyz[i];
              if(m==1) in[f] >> xyz[j];
              if(m==2) in[f] >> xyz[k];
            }
        getline(in[f], line);         /* finish the line */
      }
      else {
        getline(in[f], line_com); /* *-COORDINATE */
      }
    }

    /* << */
    out << line_com << std::endl; 
    int count=0;
    for(int k=gks-1; k<=gke-nodal; k++) 
      for(int j=gjs-1; j<=gje-nodal; j++)
        for(int i=gis-1; i<=gie-nodal; i++) {
          if(m==0) out << xyz[i] << " ";
          if(m==1) out << xyz[j] << " ";
          if(m==2) out << xyz[k] << " ";
          count++;
          if(count % 8 == 0) out << std::endl;
        }
 
    /* start each coordinate in the new line */
    out << std::endl;
  }
  xyz.resize(1);

  
  /*------------------+
  |  >> << variables  | 
  +------------------*/
  real *** var;
  alloc3d( &var, gie, gje, gke );

  for(int v=0; v<n_var; v++) { 
  
    /* >> */
    for(int f=0; f<n; f++) {
      getline(in[f], line_com);         
      if(f==0) std::cout << line_com << std::endl;  
      if(ie[f]>=is[f] && je[f]>=js[f] && ke[f]>=ks[f]) {
        for(int k=ks[f]-1; k<ke[f]; k++) 
          for(int j=js[f]-1; j<je[f]; j++)
            for(int i=is[f]-1; i<ie[f]; i++) {
              in[f] >> var[i][j][k];
            }
        getline(in[f], line);         /* finish the line */
      }
    }
  
    /* << */
    out << line_com << std::endl; 
    int count=0;
    for(int k=gks-1; k<gke; k++) 
      for(int j=gjs-1; j<gje; j++)
        for(int i=gis-1; i<gie; i++) {
          out << var[i][j][k] << " ";
          count++;
          if(count % 8 == 0) out << std::endl;
        }
  }
 
  /*---------+
  |  bodies  |
  +---------*/
  for(int io=0; ; io++) {

    /* check if the file is finished */
    bool finished = true;
    std::vector<int> ob;
    ob.resize(n);

    for(int f=0; f<n; f++) {
      in[f] >> item; 
      in[f] >> type; 
      in[f] >> item;
  
      ob[f] = atoi(item.c_str());   /* -1 if not present */

      if( item != "OF" && ob[f]>0 ) 
        finished = false;

      getline(in[f], line);         /* finish the line */
    }

    if( !finished ) {
      std::cout << "body " << io << std::endl;
 
      std::vector<Node> node_array;     /* array of nodes */
      std::vector<int>  nnodes_dom;     /* nodes in domain */
      std::vector<int>  ncells_dom;     /* cells in domain */
      std::vector<int>  cell_node;      /* cells in domain */
      std::vector<int>  first_node_dom; /* first node */
      int nnodes = 0;
      int ncells = 0;

      nnodes_dom.resize(n, 0);
      ncells_dom.resize(n, 0);
      first_node_dom.resize(n, 0);
      cell_node.resize(4);

      /*--------------------------------+
      |  get number of nodes and cells  |
      +--------------------------------*/
      for(int f=0; f<n; f++) {
        if( ob[f] >= 0 ) {
          std::cout << "present at " << f << std::endl;

          getline(in[f], line_var); 
          getline(in[f], line); 

          size_t      n1, n2;
          std::string numb;
  
          n1            = line.find("N="); n1+=2; /* n1 is the beginning of the number */
          n2            = line.find(",", n1);
          numb          = line.substr(n1,n2-n1);

          nnodes_dom[f] = atoi( numb.c_str() );

          n1            = line.find("E="); n1+=2;
          n2            = line.find(",", n1);
          numb          = line.substr(n1,n2-n1);

          ncells_dom[f] = atoi( numb.c_str() );

          ncells += ncells_dom[f]; 

          getline(in[f], store_1);             /* ZONETYPE=...      */
          if(n_var>0) getline(in[f], store_2); /* PASSIVEVARLIST... */
        }
      }

      /*----------------------+
      |  >> node coordinates  |
      +----------------------*/
      for(int f=0; f<n; f++) {
        if( ob[f] >= 0 ) {
          real x,y,z;
    
          for(int l=0; l<nnodes_dom[f]; l++) { /* read all coordinates */
            nnodes++;
            in[f] >> x >> y >> z;
            node_array.push_back( Node(nnodes, x, y, z) );
          }
        }
      }

      /*--------------+
      |  body header  |
      +--------------*/
      std::cout << "nodes " << node_array.size() << std::endl;
      std::cout << "cells " << ncells            << std::endl;

      if(type=="BODY")     out << "# BODY "     << io << std::endl; 
      out << line_var << std::endl; 
      if(type=="BODY")
        out << "ZONE T=\"TRIANGLES\", NODES=" << node_array.size() 
            << ", ELEMENTS=" << ncells << ", DATAPACKING=POINT," << std::endl;
      out << store_1  << std::endl; 
      if(n_var>0) out << store_2  << std::endl; 

      /*----------------------+
      |  << node coordinates  |
      +----------------------*/
      for(int l=0; l<node_array.size(); l++) { /* write all coordinates */
        out << node_array[l] << std::endl; 
      }

      /*--------------------------------------+
      |  estimate first node for each domain  |
      +--------------------------------------*/
      first_node_dom[0] = 0;
      for(int f=1; f<n; f++)
        first_node_dom[f] = first_node_dom[f-1] + nnodes_dom[f-1];

      /*--------------+
      |  >> << cells  |
      +--------------*/
      int nn;
      if(type=="BODY")     nn=3;
      for(int f=0; f<n; f++) {
        if( ob[f] >= 0 ) {
          for(int l=0; l<ncells_dom[f]; l++) { /* read all cells */
            /* >> */
            for(int k=0; k<nn; k++) in[f] >> cell_node[k]; 
  
            /* << */
            for(int k=0; k<nn; k++) out << cell_node[k] + first_node_dom[f] << " ";
    
            out << std::endl;           
            getline(in[f], line); /* finish the line */
          }
        }
      }
    }
    else { /* bodies have finished */
      std::cout << "bodies finished " << std::endl;
      break;
    }

  } /* io - loop through bodies */
 
  /*--------------+
  |  end of file  |
  +--------------*/
  out << "# END OF TECPLOT FILE" << std::endl; 
  out.close();
}
