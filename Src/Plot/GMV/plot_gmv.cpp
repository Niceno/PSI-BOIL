#include "plot_gmv.h"

/******************************************************************************/
void PlotGMV::plot(Domain & dm, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = & dm; // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(Body & bd, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  plot_gmv_header(nam, i);
  plot_gmv_body(bd);
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const ScalarInt & sca,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalarint(*dom, sca, "A");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec, 
                   const Scalar & sca,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_scalar(*dom, scd, "D");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec,
                   const Scalar & sca, 
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_scalar(*dom, scd, "D");
  plot_gmv_scalar(*dom, sce, "E");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Vector & vec,
                   const Scalar & sca, 
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const Scalar & scf,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_vector(*dom, vec);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_scalar(*dom, scd, "D");
  plot_gmv_scalar(*dom, sce, "E");
  plot_gmv_scalar(*dom, scf, "F");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, const Scalar & scb, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, const ScalarInt & scb,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalarint(*dom, scb, "B");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, const Scalar & scb, const Scalar & scc,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, const Scalar & scb, 
                   const Scalar & scc, const Scalar & scd,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_scalar(*dom, scd, "D");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, const Scalar & scb, const Scalar & scc,
                   const Scalar & scd, const Scalar & sce,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_scalar(*dom, scd, "D");
  plot_gmv_scalar(*dom, sce, "E");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot(const Scalar & sca, const Scalar & scb, const Scalar & scc,
                   const Scalar & scd, const Scalar & sce, const Scalar & scf,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_gmv_header(nam, i);
  plot_gmv_domain(*dom);
  plot_gmv_start_scalars();
  plot_gmv_scalar(*dom, sca, "A");
  plot_gmv_scalar(*dom, scb, "B");
  plot_gmv_scalar(*dom, scc, "C");
  plot_gmv_scalar(*dom, scd, "D");
  plot_gmv_scalar(*dom, sce, "E");
  plot_gmv_scalar(*dom, scf, "F");
  plot_gmv_end_scalars();
  plot_gmv_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotGMV::plot_gmv_header(const char * nam, const int i) {

  if( boil::cart.nproc() > 1 ) {
    std::string name_f = name_file(nam, ".gmv", i, 0);
    std::string name_l = name_file(nam, ".gmv", i, boil::cart.nproc()-1);
    boil::oout << "# Plotting: " << name_f << " ... " << name_l << boil::endl;
  } else {
    std::string name_f = name_file(nam, ".gmv", i);
    boil::oout << "# Plotting: " << name_f << boil::endl;
  }

  /* open the result file */
  std::string name = name_file(nam, ".gmv", i);
  if( boil::cart.nproc() > 1 ) 
    name   = name_file(nam, ".gmv", i, boil::cart.iam());
  out.open(name.c_str());
  
  /* header */
  out << "gmvinput ascii" << boil::endl;
}
  
/******************************************************************************/
void PlotGMV::plot_gmv_domain(const Domain & dom) { // couldn't make const 
                                           // out of "*this"
  /*--------+
  |  nodes  |
  +--------*/

  if(sh==1)
    boil::oout << "# ......... in debugging mode (with buffers)" << boil::endl; 

  int boxed_nodes=0;
  int boxed_cells=0;
  if(boxed) {
    for(int c=0; c<dom.ibody().nccells(); c++) { 
      int i,j,k;
      dom.ibody().ijk(c,&i,&j,&k);
      if(box.contains(i,j,k) || !boxed) {
        boxed_cells ++;                            
        boxed_nodes += dom.ibody().nnodes(c);
      }
    }
  }

  int isc=INT_MIN+4;  int iec=INT_MAX-4;
  int jsc=INT_MIN+4;  int jec=INT_MAX-4;
  int ksc=INT_MIN+4;  int kec=INT_MAX-4;

  if(boxed) {
    isc = box.first(Comp::i());
    iec = box.last (Comp::i());
    jsc = box.first(Comp::j());
    jec = box.last (Comp::j());
    ksc = box.first(Comp::k());
    kec = box.last (Comp::k());
  }

  int x_nodes = dom.ni()-1+2*sh;
  int y_nodes = dom.nj()-1+2*sh;
  int z_nodes = dom.nk()-1+2*sh;
  int b_nodes = dom.ibody().nnodes();

  if(boxed) {
    x_nodes = boil::mini(x_nodes, iec-isc+2);
    y_nodes = boil::mini(y_nodes, jec-jsc+2);
    z_nodes = boil::mini(z_nodes, kec-ksc+2);
    b_nodes = boil::mini(b_nodes, boxed_nodes);

    if( !box.exists() ) {
      x_nodes = 0;
      y_nodes = 0;
      z_nodes = 0;
      b_nodes = 0;
    }
  }

  int x_cells = dom.ni()-2+2*sh;
  int y_cells = dom.nj()-2+2*sh;
  int z_cells = dom.nk()-2+2*sh;
  int b_cells = dom.ibody().nccells();

  if(boxed) {
    x_cells = boil::mini(x_cells, iec-isc+1);
    y_cells = boil::mini(y_cells, jec-jsc+1);
    z_cells = boil::mini(z_cells, kec-ksc+1);
    b_cells = boil::mini(b_cells, boxed_cells);
  }

  out << "nodev " << x_nodes * y_nodes * z_nodes + b_nodes << boil::endl; 

  /* domain */
  for(int i=boil::maxi(1-sh,isc); i<boil::mini(dom.ni()+sh,iec+2); i++)
    for(int j=boil::maxi(1-sh,jsc); j<boil::mini(dom.nj()+sh,jec+2); j++)
      for(int k=boil::maxi(1-sh,ksc); k<boil::mini(dom.nk()+sh,kec+2); k++)
        out << dom.xn(i) << " " << dom.yn(j) << " " << dom.zn(k) << boil::endl;

  /* immersed boundary */
  for(int c=0; c<dom.ibody().nccells(); c++) { /* was npolys */
    int i,j,k;
    dom.ibody().ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed)
      for(int i=0; i<dom.ibody().nnodes(c); i++) 
        out << dom.ibody()[c].xn(i) << " " 
            << dom.ibody()[c].yn(i) << " " 
            << dom.ibody()[c].zn(i) << boil::endl;
  }

  /*-------------------+ 
  |  cells             |
  |                    |
  |     n3--------n7   | 
  |     /|        /|   |
  |   n1--------n5 |   |
  |    | |       | |   |
  |    | |       | |   |
  |    | |       | |   |
  |    |n2-------|n6   |
  |    |/        |/    |
  |   n0--------n4     |
  |                    |
  +-------------------*/

  out << "cells " << x_cells * y_cells * z_cells + b_cells << boil::endl; 

  boil::aout << "# ......... " 
             << x_cells * y_cells * z_cells 
             << " domain cells." << boil::endl;
  boil::aout << "# ......... " 
             << b_cells
             << " cell cuts." << boil::endl;

  /* domain */
  for(int i=1; i<x_nodes; i++)
    for(int j=1; j<y_nodes; j++) 
      for(int k=1; k<z_nodes; k++) {
        const int n0 = (i-1)* y_nodes * z_nodes  
                     + (j-1)* z_nodes  + k;
        const int n1 = n0 + 1;
        const int n2 = n0 + z_nodes;
        const int n3 = n2 + 1;
        const int n4 = n0 + y_nodes * z_nodes;
        const int n5 = n1 + y_nodes * z_nodes;
        const int n6 = n2 + y_nodes * z_nodes;
        const int n7 = n3 + y_nodes * z_nodes;
        out << "hex 8" << boil::endl;
        out << n1 << " " << n5 << " " << n7 << " " << n3 << " " 
            << n0 << " " << n4 << " " << n6 << " " << n2 << boil::endl;
      }

  /* immersed body */
  int n = x_nodes * y_nodes * z_nodes; /* nodes */
  for(int c=0; c<dom.ibody().nccells(); c++) {
    int i,j,k;
    dom.ibody().ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed) {
      out << "general 1" << boil::endl;
      out << "  " << dom.ibody().nnodes(c) << boil::endl;
      for(int i=1; i<=dom.ibody().nnodes(c); i++)
        out << n+i << " ";
      out << boil::endl;
      n += dom.ibody().nnodes(c);
    }
  }

  /*-----------+
  |  material  |
  +-----------*/
  if( dom.ibody().tpolys() > 0 ) {

    out << "material 3 0" << boil::endl;
    out << "fluid"     << boil::endl;
    out << "solid"     << boil::endl;
    out << "interface" << boil::endl;

    for(int i=1-sh; i<dom.ni()-1+sh; i++)
      for(int j=1-sh; j<dom.nj()-1+sh; j++) 
        for(int k=1-sh; k<dom.nk()-1+sh; k++) {
          if(box.contains(i,j,k) || !boxed) {
            if(dom.ibody().on(i,j,k))
              out << 1 << boil::endl;
            else
              out << 2 << boil::endl;
          }
        }

    for(int c=0; c<dom.ibody().nccells(); c++) {
      int i,j,k;
      dom.ibody().ijk(c,&i,&j,&k);
      if(box.contains(i,j,k) || !boxed) 
        out << 3 << boil::endl;
    }
  }
}

/******************************************************************************/
void PlotGMV::plot_gmv_body(const Body & bod) {
                                               
  /*--------+
  |  nodes  |
  +--------*/
  out << "nodev " << bod.nnodes() << boil::endl; 
  for(int c=0; c<bod.nccells(); c++)
    for(int i=0; i<bod.nnodes(c); i++) 
      out << bod[c].xn(i) << " " 
          << bod[c].yn(i) << " " 
          << bod[c].zn(i) << boil::endl;
                                               
  /*--------+
  |  cells  |
  +--------*/
  out << "cells " << bod.nccells() << boil::endl; 
  int n=0;
  for(int c=0; c<bod.nccells(); c++) {
    out << "general 1" << boil::endl;
    out << "  " << bod.nnodes(c) << boil::endl;
    for(int i=1; i<=bod.nnodes(c); i++)
      out << n+i << " ";
    out << boil::endl;
    n += bod.nnodes(c);
  }

  /*-----------------------+
  |  plot surface normals  |
  +-----------------------*/
  out << "velocity 0" << boil::endl;
  for(int m=0; m<3; m++)
    for(int c=0; c<bod.nccells(); c++) 
      out << bod.n(m,c) << boil::endl;
}

/******************************************************************************/
void PlotGMV::plot_gmv_start_scalars() {out << "variables" << boil::endl;}

/******************************************************************************/
void PlotGMV::plot_gmv_end_scalars() {out << "endvars" << boil::endl;}

/******************************************************************************/
void PlotGMV::plot_gmv_scalar(const Domain & dom, const Scalar & sca, 
                              const char * name) {

  /* values */
  if(sca.name().length() > 0) /* use assigned name */ 
    out << sca.name() << " 0" << boil::endl;
  else                     /* use default name */
    out << name << " 0" << boil::endl;

  /* domain */
  for(int i=1-sh; i<dom.ni()-1+sh; i++)
    for(int j=1-sh; j<dom.nj()-1+sh; j++) 
      for(int k=1-sh; k<dom.nk()-1+sh; k++) 
        if(box.contains(i,j,k) || !boxed) 
          out << sca[i][j][k] << boil::endl;

  /* plot dummy (empty cell) values on immersed body */
  for(int c=0; c<dom.ibody().nccells(); c++) {
    int i,j,k;
    dom.ibody().ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed) 
      out << sca[i][j][k] << boil::endl;
  }
}

/******************************************************************************/
void PlotGMV::plot_gmv_scalarint(const Domain & dom, const ScalarInt & sca,
                                 const char * name) {

  /* values */
  if(sca.name().length() > 0) /* use assigned name */
    out << sca.name() << " 0" << boil::endl;
  else                     /* use default name */
    out << name << " 0" << boil::endl;

  /* domain */
  for(int i=1-sh; i<dom.ni()-1+sh; i++)
    for(int j=1-sh; j<dom.nj()-1+sh; j++)
      for(int k=1-sh; k<dom.nk()-1+sh; k++)
        if(box.contains(i,j,k) || !boxed)
          out << sca[i][j][k] << boil::endl;

  /* plot dummy (empty cell) values on immersed body */
  for(int c=0; c<dom.ibody().nccells(); c++) {
    int i,j,k;
    dom.ibody().ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed)
      out << sca[i][j][k] << boil::endl;
  }
}

/******************************************************************************/
void PlotGMV::plot_gmv_vector(const Domain & dom, const Vector & vec) {

  /* values */
  out << "velocity 0" << boil::endl;
  for_m(m) {
 
    /* domain */
    for(int i=1-sh; i<dom.ni()-1+sh; i++)
      for(int j=1-sh; j<dom.nj()-1+sh; j++) 
        for(int k=1-sh; k<dom.nk()-1+sh; k++) 
          if(box.contains(i,j,k) || !boxed) {
            real vel_comp=0.0;
            if(m == Comp::u()) 
             { 
              /* interpolate as you browse */
              if( dom.ibody().off(i,j,k) ) 
                vel_comp = 0.0;
              else if( dom.ibody().cut(i,j,k) ) {
                const real fm = dom.ibody().fdxe(i,j,k);
                const real fp = dom.ibody().fdxw(i,j,k);
                if( dom.ibody().on_p(i-1,j,k) && dom.ibody().on_p(i+1,j,k) )
                  vel_comp = (fm*vec[m][i][j][k] + fp*vec[m][i+1][j][k])/(fm+fp);
                else if( dom.ibody().on_p(i-1,j,k) )
                  vel_comp = vec[m][i][j][k]*fm/(fm+fp);
                else if( dom.ibody().on_p(i+1,j,k) )
                  vel_comp = vec[m][i+1][j][k]*fp/(fm+fp);
               }
              else
                vel_comp = 0.5 * (vec[m][i][j][k] + vec[m][i+1][j][k]);
             }
            if(m == Comp::v()) 
             {                        
              /* interpolate as you browse */
              if( dom.ibody().off(i,j,k) ) 
                vel_comp = 0.0;
              else if( dom.ibody().cut(i,j,k) ) {
                const real fm = dom.ibody().fdyn(i,j,k);
                const real fp = dom.ibody().fdys(i,j,k);
                if( dom.ibody().on_p(i,j-1,k) && dom.ibody().on_p(i,j+1,k) )
                  vel_comp = (fm*vec[m][i][j][k] + fp*vec[m][i][j+1][k])/(fm+fp);
                else if( dom.ibody().on_p(i,j-1,k) ) 
                  vel_comp = vec[m][i][j][k]*fm/(fm+fp);
                else if( dom.ibody().on_p(i,j+1,k) )
                  vel_comp = vec[m][i][j+1][k]*fp/(fm+fp);
               }
              else
                vel_comp = 0.5 * (vec[m][i][j][k] + vec[m][i][j+1][k]);
             }
            if(m == Comp::w()) 
             {
              /* interpolate as you browse */
              if( dom.ibody().off(i,j,k) ) 
                vel_comp = 0.0;
              else if( dom.ibody().cut(i,j,k) ) {
                const real fm = dom.ibody().fdzt(i,j,k);
                const real fp = dom.ibody().fdzb(i,j,k);
                if( dom.ibody().on_p(i,j,k-1) && dom.ibody().on_p(i,j,k+1) )
                  vel_comp = (fm*vec[m][i][j][k]+fp*vec[m][i][j][k+1])/(fm+fp);
                else if( dom.ibody().on_p(i,j,k-1) )
                  vel_comp = vec[m][i][j][k]*fm/(fm+fp);
                else if( dom.ibody().on_p(i,j,k+1) )
                  vel_comp = vec[m][i][j][k+1]*fp/(fm+fp);
               }
              else
                vel_comp = 0.5 * (vec[m][i][j][k] + vec[m][i][j][k+1]);
             }
            if( dom.ibody().nccells()>0 )
              if( dom.ibody().fV(i,j,k) < 0.55 ) vel_comp = 0.0;
            out << vel_comp << boil::endl;
          }
        
    /* plot dummy (empty cell) values on immersed body */
    for(int c=0; c<dom.ibody().nccells(); c++) {
      int i,j,k;
      dom.ibody().ijk(c,&i,&j,&k);
      if(box.contains(i,j,k) || !boxed) 
        out << 0.0 << boil::endl;
    }
  } /* m */

}

/******************************************************************************/
void PlotGMV::plot_gmv_footer() {

  /* footer */
  out << "endgmv" << boil::endl;

  /* close a file */
  out.close();
}
