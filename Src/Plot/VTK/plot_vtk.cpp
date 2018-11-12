#include "plot_vtk.h"

/******************************************************************************/
void PlotVTK::plot(Domain & dm, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = & dm; // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(Body & bod, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  /* file name extension */
  std::string name_f = name_file(nam, ".ib.vtk", i, 0);
  std::string name_l = name_file(nam, ".ib.vtk", i, boil::cart.nproc()-1);

  if( boil::cart.nproc() > 1 ) {
    boil::oout << "# Plotting: " << name_f << " ... " << name_l << boil::endl;
  } else {
    boil::oout << "# Plotting: " << name_f << boil::endl;
  }

  /*--------------+ 
  |  open a file  |
  +--------------*/
  std::string name = name_file(nam, ".ib.vtk", i);
  if( boil::cart.nproc() > 1 ) 
    name = name_file(nam, ".ib.vtk", i, boil::cart.iam());
  out.open(name.c_str());
  
  /* header */
  out << "# vtk DataFile Version 2.0"      << boil::endl;
  out << "this is a comment use it wisely" << boil::endl;
  out << "ASCII"                           << boil::endl;
  out << "DATASET POLYDATA"                << boil::endl;

  /*--------+
  |  nodes  |
  +--------*/
  out << "POINTS " << bod.nnodes() << " float" << boil::endl; 
  for(int c=0; c<bod.nccells(); c++)
    for(int i=0; i<bod.nnodes(c); i++) 
      out << bod[c].xn(i) << " " 
          << bod[c].yn(i) << " " 
          << bod[c].zn(i) << boil::endl;
                                               
  /*--------+
  |  cells  |
  +--------*/
  int n=0; 
    for(int c=0; c<bod.nccells(); c++)
      n+=(bod.nnodes(c)+1);
  out << "POLYGONS " << bod.nccells() << " " << n << boil::endl; 
  n=0;
  for(int c=0; c<bod.nccells(); c++) {
    out << bod.nnodes(c) << " ";
    for(int i=0; i<bod.nnodes(c); i++)
      out << n+i << " ";
    out << boil::endl;
    n += bod.nnodes(c);
  }

  /*-----------------+ 
  |  close the file  |
  +-----------------*/
  out.close();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const ScalarInt & sca, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalarint(*dom, sca, "a");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot_vtk_domain(const Domain & dom) { 

  /* cell centres */
  out << boil::endl;
  out << "DATASET RECTILINEAR_GRID" << boil::endl;
  out << "DIMENSIONS " << dom.ni()-2 << " " 
                       << dom.nj()-2 << " " 
                       << dom.nk()-2 << boil::endl;

  out << "X_COORDINATES " << dom.ni()-2  << " float" << boil::endl;
  for(int i=1; i<dom.ni()-1; i++) 
    out << dom.xc(i) << boil::endl;

  out << "Y_COORDINATES " << dom.nj()-2  << " float" << boil::endl;
  for(int j=1; j<dom.nj()-1; j++) 
    out << dom.yc(j) << boil::endl;

  out << "Z_COORDINATES " << dom.nk()-2  << " float" << boil::endl;
  for(int k=1; k<dom.nk()-1; k++) 
    out << dom.zc(k) << boil::endl;
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec, 
                   const Scalar & sca,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec, 
                   const Scalar & sca,
                   const Scalar & scb,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec, 
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2)
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_scalar(*dom, scd, "d");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2)
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_scalar(*dom, scd, "d");
  plot_vtk_scalar(*dom, sce, "e");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Vector & vec,
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

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2)
      << boil::endl;
  plot_vtk_vector(*dom, vec);
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_scalar(*dom, scd, "d");
  plot_vtk_scalar(*dom, sce, "e");
  plot_vtk_scalar(*dom, scf, "f");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca,
                   const Scalar & scb,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << boil::endl;
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca,
                   const ScalarInt & scb,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << boil::endl;
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2)
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalarint(*dom, scb, "b");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << boil::endl;
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << boil::endl;
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_scalar(*dom, scd, "d");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << boil::endl;
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_scalar(*dom, scd, "d");
  plot_vtk_scalar(*dom, sce, "e");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot(const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const Scalar & scf,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  plot_vtk_header(*dom, nam, i);
  plot_vtk_domain(*dom);
  out << boil::endl;
  out << "POINT_DATA " << (dom->ni()-2)*(dom->nj()-2)*(dom->nk()-2) 
      << boil::endl;
  plot_vtk_scalar(*dom, sca, "a");
  plot_vtk_scalar(*dom, scb, "b");
  plot_vtk_scalar(*dom, scc, "c");
  plot_vtk_scalar(*dom, scd, "d");
  plot_vtk_scalar(*dom, sce, "e");
  plot_vtk_scalar(*dom, scf, "f");
  plot_vtk_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotVTK::plot_vtk_header(const Domain & dom, 
                              const char * nam, const int i) {

  if( boil::cart.nproc() > 1 ) {
    std::string name_f = name_file(nam, ".vtk", i, 0);
    std::string name_l = name_file(nam, ".vtk", i, boil::cart.nproc()-1);
    boil::oout << "# Plotting: " << name_f << " ... " << name_l << boil::endl;
  } else {
    std::string name_f = name_file(nam, ".vtk", i);
    boil::oout << "# Plotting: " << name_f << boil::endl;
  }


  /* open a file */
  std::string name = name_file(nam, ".vtk", i);
  if( boil::cart.nproc() > 1 ) 
    name = name_file(nam, ".vtk", i, boil::cart.iam());
  out.open(name.c_str());
  
  /* header */
  out << "# vtk DataFile Version 2.0" << boil::endl;

  Range<int> c_rang_i = dom.cxg();
  Range<int> c_rang_j = dom.cyg();
  Range<int> c_rang_k = dom.czg();

  out << "RANGES " << c_rang_i.first() << " " 
                   << c_rang_i.last()  << " "
                   << c_rang_j.first() << " " 
                   << c_rang_j.last()  << " "
                   << c_rang_k.first() << " " 
                   << c_rang_k.last()  << boil::endl;

  out << "ASCII" << boil::endl;
}
  
/******************************************************************************/
void PlotVTK::plot_vtk_scalar(const Domain & dom, const Scalar & sca, 
                              const char * nm) {

  /* values */
  out << boil::endl;

  if(sca.name().length() > 0) /* use assigned name */
    out << "SCALARS " << sca.name() << " float" << boil::endl;
  else                     /* use default name */
    out << "SCALARS " << nm << " float" << boil::endl;

  /* values */
  out << "LOOKUP_TABLE default" << boil::endl;

  for(int k=1; k<dom.nk()-1; k++) 
    for(int j=1; j<dom.nj()-1; j++) 
      for(int i=1; i<dom.ni()-1; i++)
        out << sca[i][j][k] << boil::endl;
}

/******************************************************************************/
void PlotVTK::plot_vtk_scalarint(const Domain & dom, const ScalarInt & sca,
                                 const char * nm) {

  /* values */
  out << boil::endl;

  if(sca.name().length() > 0) /* use assigned name */
    out << "SCALARS " << sca.name() << " float" << boil::endl;
  else                     /* use default name */
    out << "SCALARS " << nm << " float" << boil::endl;

  /* values */
  out << "LOOKUP_TABLE default" << boil::endl;

  for(int k=1; k<dom.nk()-1; k++)
    for(int j=1; j<dom.nj()-1; j++)
      for(int i=1; i<dom.ni()-1; i++)
        out << sca[i][j][k] << boil::endl;
}

/******************************************************************************/
void PlotVTK::plot_vtk_vector(const Domain & dom, const Vector & vec) {

  /* values */
  out << boil::endl;
  out << "VECTORS velocity float" << boil::endl;

  Comp m;
  for(int k=1; k<dom.nk()-1; k++) 
    for(int j=1; j<dom.nj()-1; j++) 
      for(int i=1; i<dom.ni()-1; i++) {
        m=Comp::u(); real u = 0.5 * (vec[m][i][j][k] + vec[m][i+1][j][k]);
        m=Comp::v(); real v = 0.5 * (vec[m][i][j][k] + vec[m][i][j+1][k]);
        m=Comp::w(); real w = 0.5 * (vec[m][i][j][k] + vec[m][i][j][k+1]);

        out << u << " " << v << " " << w << boil::endl;
      }	
}

/******************************************************************************/
void PlotVTK::plot_vtk_footer() {

  /* close a file */
  out.close();
}
