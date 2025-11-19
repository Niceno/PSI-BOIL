#include "connect.h"

/******************************************************************************/
void connect_vtk(const char * bname, const int n, const int t) {

  /*--------------------+
  |  working variables  |
  +--------------------*/
  std::string item;
  std::string name;
  std::string line;
  std::vector<real> x;
  std::vector<real> y;
  std::vector<real> z;

  /*---------------+
  |  file streams  |
  +---------------*/
  std::ofstream out;
  std::ifstream * in = new std::ifstream[n];

  /*---------------------------------------+
  |  >> browse through files to open them  |
  +---------------------------------------*/
  for(int f=0; f<n; f++) {
    std::string current = name_file(bname, ".vtk", t, f);
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
  int gis=2, gie=0, gjs=2, gje=0, gks=2, gke=0;
  int nodal=0;

  for(int f=0; f<n; f++) {
    getline(in[f], line);         /* "# vtk DataFile Version 2.0" */

    /* get ranges */
    in[f] >> item; in[f] >> is[f]; in[f] >> ie[f];
                   in[f] >> js[f]; in[f] >> je[f];
                   in[f] >> ks[f]; in[f] >> ke[f];

    if(is[f] < gis) gis=is[f]; if(ie[f] > gie) gie=ie[f];
    if(js[f] < gjs) gjs=js[f]; if(je[f] > gje) gje=je[f];
    if(ks[f] < gks) gks=ks[f]; if(ke[f] > gke) gke=ke[f];

    getline(in[f], line);         /* finish the line */
  }
  
  /*---------------------------+
  |  resize coordinate arrays  |
  +---------------------------*/
  x.resize(gie+1);
  y.resize(gje+1);
  z.resize(gke+1);

  /*---------------+
  |  >> grid info  | 
  +---------------*/
  int n_var;
  for(int f=0; f<n; f++) {
    getline(in[f], line);     /* "ASCII" */
    getline(in[f], line);     /* "" */
    getline(in[f], line);     /* "DATASET RECTILINEAR_GRID" */
    getline(in[f], line);     /* "DIMENSIONS ..." */
    /*----+
    |  x  |
    +----*/
    getline(in[f], line);     /* "X_COORDINATES ..." */
    for(int i=is[f]; i<=ie[f]; i++) {
      in[f] >> x[i];     
    }
    getline(in[f], line);     /* finish the line */
    /*----+
    |  y  |
    +----*/
    getline(in[f], line);     /* "Y_COORDINATES ..." */
    for(int j=js[f]; j<=je[f]; j++) {
      in[f] >> y[j];     
    }
    getline(in[f], line);     /* finish the line */
    /*----+
    |  z  |
    +----*/
    getline(in[f], line);     /* "Z_COORDINATES ..." */
    for(int k=ks[f]; k<=ke[f]; k++) {
      in[f] >> z[k];     
    }
    getline(in[f], line);     /* finish the line */
  }

  /*----------------------+
  |  << start the output  |
  +----------------------*/
  std::string out_name = name_file(bname, ".vtk", t);
  std::cout << "creating: " << out_name << std::endl;
  out.open( out_name.c_str() );

  out << "# vtk DataFile Version 2.0"                       << std::endl;
  out << "RANGES 1 " << gie << " 1 " << gje << " 1 " << gke << std::endl;
  out << "ASCII"                                            << std::endl;
  out << std::endl;
  out << "DATASET RECTILINEAR_GRID"                         << std::endl;
  out << "DIMENSIONS " << gie << " " << gje << " " << gke   << std::endl;

  /*----------------------+
  |  << node coordinates  | 
  +----------------------*/
  out << "X_COORDINATES " << gie << " float" << std::endl;
  for(int i=1; i<=gie; i++) out << x[i]      << std::endl;
  out << "Y_COORDINATES " << gje << " float" << std::endl;
  for(int j=1; j<=gje; j++) out << y[j]      << std::endl;
  out << "Z_COORDINATES " << gke << " float" << std::endl;
  for(int k=1; k<=gke; k++) out << z[k]      << std::endl;

  /*----------------+
  |  >> point data  |
  +----------------*/
  for(int f=0; f<n; f++) {
    in[f] >> item;
    getline(in[f], line);     /* "POINT DATA ..." */
  }

  /*------------------+
  |  >> << variables  | 
  +------------------*/
  if( item=="POINT_DATA" ) {
    real *** phi;
    real *** u;
    real *** v;
    real *** w;
    alloc3d( &phi, gie, gje, gke );
    alloc3d( &u,   gie, gje, gke );
    alloc3d( &v,   gie, gje, gke );
    alloc3d( &w,   gie, gje, gke );

    out << "POINT_DATA " << gie * gje * gke << std::endl;

    for(int cnt=0; cnt<10; cnt++) {

      /*-----+
      |  >>  |
      +-----*/
      for(int f=0; f<n; f++) {
   
        item = "END";

        getline(in[f], line);     /* "" */
        in[f] >> item;
        in[f] >> name;
        if(f==0 && item!="END") std::cout << "processing " << name << std::endl;
        getline(in[f], line);     /* finish the line */

        if( item=="VECTORS" ) {
          for(int k=ks[f]-1; k<ke[f]; k++) 
            for(int j=js[f]-1; j<je[f]; j++)
              for(int i=is[f]-1; i<ie[f]; i++) {
                in[f] >> u[i][j][k];
                in[f] >> v[i][j][k];
                in[f] >> w[i][j][k];
              }
        }
        else if( item=="SCALARS" ) {
          getline(in[f], line);     /* "LOOKUP_TABLE default" */
          for(int k=ks[f]-1; k<ke[f]; k++) 
            for(int j=js[f]-1; j<je[f]; j++)
              for(int i=is[f]-1; i<ie[f]; i++) {
                in[f] >> phi[i][j][k];
              }
          getline(in[f], line);     /* finish the line */
        }
        else {
          /*--------------+
          |  end of file  |
          +--------------*/
          out.close();
          exit(0);
        } 
      }

      /*-----+
      |  <<  |
      +-----*/
      if( item=="VECTORS" ) {
        out << std::endl;
        out << "VECTORS vec float" << std::endl;
        for(int k=0; k<gke; k++) 
          for(int j=0; j<gje; j++) 
            for(int i=0; i<gie; i++) 
              out << u[i][j][k] << " " 
                  << v[i][j][k] << " " 
                  << w[i][j][k] << std::endl;
      }
      if( item=="SCALARS" ) {
        out << std::endl;
        out << "SCALARS " << name << " float" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;
        for(int k=0; k<gke; k++) 
          for(int j=0; j<gje; j++) 
            for(int i=0; i<gie; i++) 
              out << phi[i][j][k] << std::endl;
      }

    }
  }
}
