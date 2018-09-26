#include "plot_tec.h"

/******************************************************************************/
void PlotTEC::plot(Domain & dm, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = & dm; // take it as a constant

  std::vector<int> vars;

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(Body & bd, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i) {

  std::vector<int> vars;

  boil::timer.start("plotting");
  plot_tec_header(nam, i);
  plot_tec_body(bd, vars);
  plot_tec_footer();
  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; vars.push_back(4);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("PHI"); if(sca.name().length() > 0) vnames[3] = sca.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const ScalarInt & sca, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; vars.push_back(4);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("PHI"); if(sca.name().length() > 0) vnames[3] = sca.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalarint(*dom, sca);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=6; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec, 
                   const Scalar & sca,
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=7; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames [6] = sca.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=8; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames [6] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames [7] = scb.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=9; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames [6] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames [7] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames [8] = scc.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=10; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames [6] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames [7] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames [8] = scc.name();
  vnames.push_back("D"); if(scd.name().length() > 0) vnames [9] = scd.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=11; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames [6] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames [7] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames [8] = scc.name();
  vnames.push_back("D"); if(scd.name().length() > 0) vnames [9] = scd.name();
  vnames.push_back("E"); if(sce.name().length() > 0) vnames[10] = sce.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec,
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

  std::vector<int> vars;
  for(int v=4; v<=12; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("U"); 
  vnames.push_back("V"); 
  vnames.push_back("W"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames [6] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames [7] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames [8] = scc.name();
  vnames.push_back("D"); if(scd.name().length() > 0) vnames [9] = scd.name();
  vnames.push_back("E"); if(sce.name().length() > 0) vnames[10] = sce.name();
  vnames.push_back("F"); if(scf.name().length() > 0) vnames[11] = scf.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_scalar(*dom, scf);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const Scalar & scb, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=5; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames[3] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames[4] = scb.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca,
                   const ScalarInt & scb,
                   const char * nam,
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=5; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X");
  vnames.push_back("Y");
  vnames.push_back("Z");
  vnames.push_back("A"); if(sca.name().length() > 0) vnames[3] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames[4] = scb.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalarint(*dom, scb);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const Scalar & scb, 
                   const Scalar & scc, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=6; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames[3] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames[4] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames[5] = scc.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const Scalar & scb, 
                   const Scalar & scc, 
                   const Scalar & scd, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=7; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames[3] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames[4] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames[5] = scc.name();
  vnames.push_back("D"); if(scd.name().length() > 0) vnames[6] = scd.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const Scalar & scb, 
                   const Scalar & scc, 
                   const Scalar & scd, 
                   const Scalar & sce, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=8; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames[3] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames[4] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames[5] = scc.name();
  vnames.push_back("D"); if(scd.name().length() > 0) vnames[6] = scd.name();
  vnames.push_back("E"); if(sce.name().length() > 0) vnames[7] = sce.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const Scalar & scb, 
                   const Scalar & scc, 
                   const Scalar & scd, 
                   const Scalar & sce, 
                   const Scalar & scf, 
                   const char * nam, 
                   const int i) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; 
  for(int v=4; v<=9; v++) vars.push_back(v);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("A"); if(sca.name().length() > 0) vnames[3] = sca.name();
  vnames.push_back("B"); if(scb.name().length() > 0) vnames[4] = scb.name();
  vnames.push_back("C"); if(scc.name().length() > 0) vnames[5] = scc.name();
  vnames.push_back("D"); if(scd.name().length() > 0) vnames[6] = scd.name();
  vnames.push_back("E"); if(sce.name().length() > 0) vnames[7] = sce.name();
  vnames.push_back("F"); if(scf.name().length() > 0) vnames[8] = scf.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_scalar(*dom, scf);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot_tec_header(const Domain & dom, 
                              const char * nam, const int i) {

  if( boil::cart.nproc() > 1 ) {
    std::string name_f = name_file(nam, ".dat", i, 0);
    std::string name_l = name_file(nam, ".dat", i, boil::cart.nproc()-1);
    boil::oout << "# Plotting: " << name_f << " ... " << name_l << boil::endl;
  } else {
    std::string name_f = name_file(nam, ".dat", i);
    boil::oout << "# Plotting: " << name_f << boil::endl;
  }

  if(sh==1)
    boil::oout << "# ......... in debugging mode (with buffers)" << boil::endl; 

  /* open the result file */
  std::string name = name_file(nam, ".dat", i);
  if( boil::cart.nproc() > 1 ) 
    name = name_file(nam, ".dat", i, boil::cart.iam());
  out.open(name.c_str());
  out << "# FILE CREATED WITH PSI-BOIL" << boil::endl;

  if(nodal) out << "# NODAL FORMAT " << boil::endl;
  else 
    if(sh==1){
      out << "# CELLCENTERED FORMAT WITH BUFFER" << boil::endl;
    } else {
      out << "# CELLCENTERED FORMAT " << boil::endl;
    }


  Range<int> c_rang_i = dom.cxg();
  Range<int> c_rang_j = dom.cyg();
  Range<int> c_rang_k = dom.czg();

  if(boxed) {
    c_rang_i.first( box.first(Comp::i()) );
    c_rang_i.last ( box.last (Comp::i()) );
    c_rang_j.first( box.first(Comp::j()) );
    c_rang_j.last ( box.last (Comp::j()) );
    c_rang_k.first( box.first(Comp::k()) );
    c_rang_k.last ( box.last (Comp::k()) );
  }

  if(boxed) {
    out << "# I-RANGE " << dom.global_I(c_rang_i.first()) << " " 
                        << dom.global_I(c_rang_i.last())  << boil::endl;
    out << "# J-RANGE " << dom.global_J(c_rang_j.first()) << " " 
                        << dom.global_J(c_rang_j.last())  << boil::endl;
    out << "# K-RANGE " << dom.global_K(c_rang_k.first()) << " " 
                        << dom.global_K(c_rang_k.last())  << boil::endl;
  } else {
    out << "# I-RANGE " << c_rang_i.first() << " " 
                        << c_rang_i.last()  << boil::endl;
    out << "# J-RANGE " << c_rang_j.first() << " " 
                        << c_rang_j.last()  << boil::endl;
    out << "# K-RANGE " << c_rang_k.first() << " " 
                        << c_rang_k.last()  << boil::endl;
  }
}
  
/******************************************************************************/
void PlotTEC::plot_tec_header(const char * nam, const int i) {

  /* file name extension */
  std::string name_f = name_file(nam, ".dat", i, 0);
  std::string name_l = name_file(nam, ".dat", i, boil::cart.nproc()-1);

  if( boil::cart.nproc() > 1 ) {
    boil::oout << "# Plotting: " << name_f << " ... " << name_l << boil::endl;
  } else {
    boil::oout << "# Plotting: " << name_f << boil::endl;
  }

  if(sh==1)
    boil::oout << "# ......... in debugging mode (with buffers)" << boil::endl; 

  /* open the result file */
  std::string name = name_file(nam, ".dat", i);
  if( boil::cart.nproc() > 1 ) 
    name = name_file(nam, ".dat", i, boil::cart.iam());
  out.open(name.c_str());

  out << "# FILE CREATED WITH PSI-BOIL" << boil::endl;

}
  
/******************************************************************************/
void PlotTEC::plot_tec_prologue(const std::vector<std::string> & vnames) {

  int x_nodes = dom->ni()-1+2*sh;
  int y_nodes = dom->nj()-1+2*sh;
  int z_nodes = dom->nk()-1+2*sh;

  if(boxed) {
    int isc = box.first(Comp::i());
    int iec = box.last (Comp::i());
    int jsc = box.first(Comp::j());
    int jec = box.last (Comp::j());
    int ksc = box.first(Comp::k());
    int kec = box.last (Comp::k());

    if( !box.exists() ) {
      x_nodes = 0;
      y_nodes = 0;
      z_nodes = 0;
    }

    x_nodes = boil::mini(x_nodes, iec-isc+2);
    y_nodes = boil::mini(y_nodes, jec-jsc+2);
    z_nodes = boil::mini(z_nodes, kec-ksc+2);
  }

  out << "VARIABLES=";         
  for(int i=0; i<vnames.size(); i++) out << "\"" << vnames[i] << "\" " ;
  out << boil::endl;
  out << " ZONE I=" << x_nodes - nodal
      <<     ", J=" << y_nodes - nodal
      <<     ", K=" << z_nodes - nodal;
  out << " DATAPACKING=BLOCK";
  if(vnames.size() == 4) {
    if(nodal) out << ", VARLOCATION=([4]=NODAL)";
    else      out << ", VARLOCATION=([4]=CELLCENTERED)";
  }
  if(vnames.size() > 4) {
    if(nodal) out << ", VARLOCATION=([4-" << vnames.size() << "]=NODAL)";
    else      out << ", VARLOCATION=([4-" << vnames.size() << "]=CELLCENTERED)";
  }
  out << boil::endl;
}

/******************************************************************************/
void PlotTEC::plot_tec_domain(const Domain & dom) {

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

  for_m(m) {
    if(m==Comp::u()) out << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) out << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) out << "# Z-COORDINATES" << boil::endl;
    int count=0;
    if(nodal) {
      for(int k=boil::maxi(1-sh,ksc); 
              k<boil::mini(dom.nk()+sh-nodal,kec+2-nodal); k++) 
        for(int j=boil::maxi(1-sh,jsc); 
                j<boil::mini(dom.nj()+sh-nodal,jec+2-nodal); j++)
          for(int i=boil::maxi(1-sh,isc); 
                  i<boil::mini(dom.ni()+sh-nodal,iec+2-nodal); i++) {
            if(m==Comp::u()) out << 0.5*(dom.xn(i)+dom.xn(i+1)) << " ";
            if(m==Comp::v()) out << 0.5*(dom.yn(j)+dom.yn(j+1)) << " ";
            if(m==Comp::w()) out << 0.5*(dom.zn(k)+dom.zn(k+1)) << " ";

            count++;
            if(count % 8 == 0) out << boil::endl;
          }
    } else {
      for(int k=boil::maxi(1-sh,ksc); 
              k<boil::mini(dom.nk()+sh,kec+2); k++) 
        for(int j=boil::maxi(1-sh,jsc); 
                j<boil::mini(dom.nj()+sh,jec+2); j++)
          for(int i=boil::maxi(1-sh,isc); 
                  i<boil::mini(dom.ni()+sh,iec+2); i++) {
            if(m==Comp::u()) out << dom.xn(i) << " ";
            if(m==Comp::v()) out << dom.yn(j) << " ";
            if(m==Comp::w()) out << dom.zn(k) << " ";

            count++;
            if(count % 8 == 0) out << boil::endl;
          }
    }
 
     /* start each coordinate in the new line (if needed) */
     if(count % 8 != 0) out << boil::endl;
  } /* m */
}

/******************************************************************************/
void PlotTEC::plot_tec_body(const Body & bod, std::vector<int> & vars) {

  if( bod.tpolys() < 1) 
    return;

  /* body's numbers of nodes and cells */
  int b_nodes = bod.nnodes();
  int b_cells = bod.nccells();

  /* correct b_nodes and b_cells if boxed */
  if(boxed) {
    int boxed_nodes=0;
    int boxed_cells=0;
    for(int c=0; c<bod.nccells(); c++) { 
      int i,j,k;
      bod.ijk(c,&i,&j,&k);
      if(box.contains(i,j,k)) {
        boxed_cells ++;                            
        boxed_nodes += bod.nnodes(c);
      }
    }
    b_nodes = boil::mini(b_nodes, boxed_nodes);
    b_cells = boil::mini(b_cells, boxed_cells);
  }

  if(b_cells < 1) {
    out << "# BODY " << -1 << boil::endl;
    return;
  }

  /*-------------------------------------------+
  |  count the number of additional triangles  |
  +-------------------------------------------*/
  int nat=0; /* number of additional triangles */
  for(int c=0; c<bod.nccells(); c++) {
    assert( bod.nnodes(c) >= 3 );
    int i,j,k;
    bod.ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed) {
      nat += bod.nnodes(c)-3;
    }
  } 

  /*----------------------+
  |  header for the zone  |
  +----------------------*/
  out << "# BODY " << 1 << boil::endl;
  out << "VARIABLES=\"X\" \"Y\" \"Z\" " << boil::endl;
  out << "ZONE T=\"TRIANGLES\", N=" << b_nodes
      << ", E=" << b_cells + nat       
      << ", DATAPACKING=POINT," << boil::endl;
  out << "ZONETYPE=FETRIANGLE" << boil::endl;
  if(vars.size() > 0) {
    out << "PASSIVEVARLIST=[";
    for(int i=0; i<vars.size(); i++) {
      if(i==0) 
        out << vars[i];
      else
        out << "," << vars[i];
    }
    out << "]" << boil::endl;
  }

  /*--------+
  |  nodes  |
  +--------*/
  for(int c=0; c<bod.nccells(); c++) {
    int i,j,k;
    bod.ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed) 
      for(int i=0; i<bod.nnodes(c); i++) {
        out << bod[c].xn(i) << " " 
            << bod[c].yn(i) << " " 
            << bod[c].zn(i) << std::endl;
    }
  }
                                               
  /*--------+
  |  cells  |
  +--------*/
  int n=0;
  for(int c=0; c<bod.nccells(); c++) {
    int i,j,k;
    bod.ijk(c,&i,&j,&k);
    if(box.contains(i,j,k) || !boxed) {
      int nn = bod.nnodes(c);
                out << n+1 << " " << n+2 << " " << n+3 << boil::endl;
      if(nn>3)  out << n+1 << " " << n+3 << " " << n+4 << boil::endl;
      if(nn>4)  out << n+1 << " " << n+4 << " " << n+5 << boil::endl;
      if(nn>5)  out << n+1 << " " << n+5 << " " << n+6 << boil::endl;
      n+=nn;
    }
  }
}

/******************************************************************************/
void PlotTEC::plot_tec_scalar(const Domain & dom, const Scalar & sca) {

  out << "# SCALAR VALUES" << boil::endl;

  int count=0;
  for(int k=1-sh; k<dom.nk()-1+sh; k++) 
    for(int j=1-sh; j<dom.nj()-1+sh; j++) 
      for(int i=1-sh; i<dom.ni()-1+sh; i++) {
        if(box.contains(i,j,k) || !boxed) {
          out << sca[i][j][k] << " ";       
          count++;
          if(count % 8 == 0) out << boil::endl;
        }
      }
 
  /* end with new line (if needed) */
  if(count % 8 != 0) out << boil::endl;
}

/******************************************************************************/
void PlotTEC::plot_tec_scalarint(const Domain & dom, const ScalarInt & sca) {

  out << "# SCALAR VALUES" << boil::endl;

  int count=0;
  for(int k=1-sh; k<dom.nk()-1+sh; k++) 
    for(int j=1-sh; j<dom.nj()-1+sh; j++) 
      for(int i=1-sh; i<dom.ni()-1+sh; i++) {
        if(box.contains(i,j,k) || !boxed) {
          out << sca[i][j][k] << " ";       
          count++;
          if(count % 8 == 0) out << boil::endl;
        }
      }
 
  /* end with new line (if needed) */
  if(count % 8 != 0) out << boil::endl;
}

/******************************************************************************/
void PlotTEC::plot_tec_vector(const Domain & dom, const Vector & vec) {

  for_m(m) {
    if(m==Comp::u()) out << "# U-VELOCITY" << boil::endl;
    if(m==Comp::v()) out << "# V-VELOCITY" << boil::endl;
    if(m==Comp::w()) out << "# W-VELOCITY" << boil::endl;
    int count=0;
    for(int k=1-sh; k<dom.nk()-1+sh; k++) 
      for(int j=1-sh; j<dom.nj()-1+sh; j++)
        for(int i=1-sh; i<dom.ni()-1+sh; i++) {
          if(box.contains(i,j,k) || !boxed) {
            real vel_comp;
            if(m == Comp::u()) 
             { 
              /* interpolate as you browse */
              if( dom.ibody().off(i,j,k) ) 
                vel_comp = 0.0;
              else if( dom.ibody().cut(i,j,k) ) {
                const real fm = boil::minr(0.5*dom.dxc(i)
                                          ,dom.ibody().fdxe(i,j,k)*dom.dxe(i));
                const real fp = boil::minr(0.5*dom.dxc(i)
                                          ,dom.ibody().fdxw(i,j,k)*dom.dxw(i));
                vel_comp = (vec[m][i][j][k]*fm+vec[m][i+1][j][k]*fp) / (fm+fp);
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
                const real fm = boil::minr(0.5*dom.dyc(j)
                                          ,dom.ibody().fdyn(i,j,k)*dom.dyn(j));
                const real fp = boil::minr(0.5*dom.dyc(j)
                                          ,dom.ibody().fdys(i,j,k)*dom.dys(j));
                vel_comp = (vec[m][i][j][k]*fm+vec[m][i][j+1][k]*fp) / (fm+fp);
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
                const real fm = boil::minr(0.5*dom.dzc(k)
                                          ,dom.ibody().fdzt(i,j,k)*dom.dzt(k));
                const real fp = boil::minr(0.5*dom.dzc(k)
				          ,dom.ibody().fdzb(i,j,k)*dom.dzb(k));
                vel_comp = (vec[m][i][j][k]*fm+vec[m][i][j][k+1]*fp) / (fm+fp);
               }
              else
                vel_comp = 0.5 * (vec[m][i][j][k] + vec[m][i][j][k+1]);
             }
            out << vel_comp << " ";
            count++;
            if(count % 8 == 0) out << boil::endl;
          } /* boxed */
        }
 
    /* start each coordinate in the new line (if needed) */
    if(count % 8 != 0) out << boil::endl;
   }
}

/******************************************************************************/
void PlotTEC::plot_tec_footer() {

  out << "# END OF TECPLOT FILE" << boil::endl;

  /* close the file */
  out.close();
}

/*-----------------------------------------------------------------------------+
 '$Id: plot_tec.cpp,v 1.44 2015/05/05 14:27:20 sato Exp $'/
+-----------------------------------------------------------------------------*/
