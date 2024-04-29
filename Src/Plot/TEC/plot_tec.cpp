#include "plot_tec.h"

/******************************************************************************/
void PlotTEC::plot(Domain & dm, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = & dm; // take it as a constant

  std::vector<int> vars;

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames, t);
  plot_tec_domain(*dom);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(Body & bd, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i, Times * t) {

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
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; vars.push_back(4);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("PHI"); if(sca.name().length() > 0) vnames[3] = sca.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames, t);
  plot_tec_domain(*dom);
  plot_tec_scalar(*dom, sca);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const ScalarInt & sca, 
                   const char * nam, 
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = sca.domain(); // take it as a constant

  std::vector<int> vars; vars.push_back(4);

  std::vector<std::string> vnames;
  vnames.push_back("X"); 
  vnames.push_back("Y"); 
  vnames.push_back("Z"); 
  vnames.push_back("PHI"); if(sca.name().length() > 0) vnames[3] = sca.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames, t);
  plot_tec_domain(*dom);
  plot_tec_scalarint(*dom, sca);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Vector & vec, 
                   const char * nam, 
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i,
                   Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i,
                   Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
void PlotTEC::plot(const Vector & vec,
                   const Scalar & sca,
                   const Scalar & scb,
                   const Scalar & scc,
                   const Scalar & scd,
                   const Scalar & sce,
                   const Scalar & scf,
                   const Scalar & scg,
                   const char * nam,
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=13; v++) vars.push_back(v);

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
  vnames.push_back("G"); if(scg.name().length() > 0) vnames[12] = scg.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames, t);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_scalar(*dom, scf);
  plot_tec_scalar(*dom, scg);
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
                   const Scalar & scg,
                   const Scalar & sch,
                   const char * nam,
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=14; v++) vars.push_back(v);

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
  vnames.push_back("G"); if(scg.name().length() > 0) vnames[12] = scg.name();
  vnames.push_back("H"); if(sch.name().length() > 0) vnames[13] = sch.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames, t);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_scalar(*dom, scf);
  plot_tec_scalar(*dom, scg);
  plot_tec_scalar(*dom, sch);
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
                   const Scalar & scg,
                   const Scalar & sch,
                   const Scalar & sci,
                   const char * nam,
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = vec.domain(); // take it as a constant

  std::vector<int> vars;
  for(int v=4; v<=15; v++) vars.push_back(v);

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
  vnames.push_back("G"); if(scg.name().length() > 0) vnames[12] = scg.name();
  vnames.push_back("H"); if(sch.name().length() > 0) vnames[13] = sch.name();
  vnames.push_back("I"); if(sci.name().length() > 0) vnames[14] = sci.name();

  plot_tec_header(*dom, nam, i);
  plot_tec_prologue(vnames, t);
  plot_tec_domain(*dom);
  plot_tec_vector(*dom, vec);
  plot_tec_scalar(*dom, sca);
  plot_tec_scalar(*dom, scb);
  plot_tec_scalar(*dom, scc);
  plot_tec_scalar(*dom, scd);
  plot_tec_scalar(*dom, sce);
  plot_tec_scalar(*dom, scf);
  plot_tec_scalar(*dom, scg);
  plot_tec_scalar(*dom, sch);
  plot_tec_scalar(*dom, sci);
  plot_tec_body(dom->ibody(),vars);
  plot_tec_footer();

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTEC::plot(const Scalar & sca, 
                   const Scalar & scb, 
                   const char * nam, 
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
                   const int i, Times * t) {

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
  plot_tec_prologue(vnames, t);
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
void PlotTEC::plot(const Pathline & pl, 
                   const char * nam, 
                   const int i, Times * t) {
  if( boil::cart.iam()==0 ) {
    /* open the result file */
    std::string name = name_file(nam, ".dat", i);
    boil::oout << "# Plotting: "<<name<<"\n";
    out.open(name.c_str());

    // visit cannot read the next line.
    //out << "#Number_of_variables= "<<6+pl.nval()<<"\n";
    out << "VARIABLES= X Y Z U V W ";
    std::string vname;
    if (pl.nval()>=1) { 
      vname = "A"; if(pl.s1->name().length() > 0) vname = pl.s1->name();
      out << "\""<< vname <<"\" " ;
    }
    if (pl.nval()>=2) { 
      vname = "B"; if(pl.s2->name().length() > 0) vname = pl.s2->name();
      out << "\""<< vname <<"\" " ;
    }
    if (pl.nval()>=3) { 
      vname = "C"; if(pl.s3->name().length() > 0) vname = pl.s3->name();
      out << "\""<< vname <<"\" " ;
    }
    out << boil::endl;

    out << "ZONE I= " <<pl.np()<<" DATAPACKING=POINT\n";
    if (t!=NULL) {
      out << "SOLUTIONTIME= "<<t->current_time()<<"\n";
    }
    for (int ip = 0; ip < pl.np(); ip++){
      out << pl.particles[ip].x()<<" "
          << pl.particles[ip].y()<<" "
          << pl.particles[ip].z()<<" "
          << pl.particles[ip].u()<<" "
          << pl.particles[ip].v()<<" "
          << pl.particles[ip].w()<<" ";
      //if (pl.nval()>=1) out << pl.particles[ip].sval(1)<<" ";
      //if (pl.nval()>=2) out << pl.particles[ip].s2()<<" ";
      //if (pl.nval()>=3) out << pl.particles[ip]->sval(3)<<" ";
      for (int ival = 0; ival < pl.nval(); ival++) {
        out << pl.particles[ip].sval(ival)<<" ";
      }
      out << boil::endl;
    }
    out.close();
  }
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

  out << "# I-RANGE " << c_rang_i.first() << " " 
                      << c_rang_i.last()  << boil::endl;
  out << "# J-RANGE " << c_rang_j.first() << " " 
                      << c_rang_j.last()  << boil::endl;
  out << "# K-RANGE " << c_rang_k.first() << " " 
                      << c_rang_k.last()  << boil::endl;
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
void PlotTEC::plot_tec_prologue(const std::vector<std::string> & vnames,
                                const Times * t) {

  int x_nodes = dom->ni()-2*boil::BW+1 + 2*sh;
  int y_nodes = dom->nj()-2*boil::BW+1 + 2*sh;
  int z_nodes = dom->nk()-2*boil::BW+1 + 2*sh;

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
  if(t != NULL){
    out << "SOLUTIONTIME= "<<t->current_time()<<"\n";
  }
}

/******************************************************************************/
void PlotTEC::plot_tec_domain(const Domain & dom) {

  for_m(m) {
    if(m==Comp::u()) out << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) out << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) out << "# Z-COORDINATES" << boil::endl;
    int count=0;
    if(nodal) {
      for(int k=boil::BW-sh; k<dom.nk()-boil::BW+sh-nodal+1; k++) 
        for(int j=boil::BW-sh; j<dom.nj()-boil::BW+sh-nodal+1; j++)
          for(int i=boil::BW-sh; i<dom.ni()-boil::BW+sh-nodal+1; i++) {
            if(m==Comp::u()) out << 0.5*(dom.xn(i)+dom.xn(i+1)) << " ";
            if(m==Comp::v()) out << 0.5*(dom.yn(j)+dom.yn(j+1)) << " ";
            if(m==Comp::w()) out << 0.5*(dom.zn(k)+dom.zn(k+1)) << " ";

            count++;
            if(count % 8 == 0) out << boil::endl;
          }
    } else {
      for(int k=boil::BW-sh; k<dom.nk()-boil::BW+sh+1; k++) 
        for(int j=boil::BW-sh; j<dom.nj()-boil::BW+sh+1; j++)
          for(int i=boil::BW-sh; i<dom.ni()-boil::BW+sh+1; i++) {
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

  if(b_cells < 1) {
    out << "# BODY " << -1 << boil::endl;
    return;
  }

  if(!b_plot_body) return;

  /*-------------------------------------------+
  |  count the number of additional triangles  |
  +-------------------------------------------*/
  int nat=0; /* number of additional triangles */
  for(int c=0; c<bod.nccells(); c++) {
    assert( bod.nnodes(c) >= 3 );
    int i,j,k;
    bod.ijk(c,&i,&j,&k);
    nat += bod.nnodes(c)-3;
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
    //if(box.contains(i,j,k) || !boxed)
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
    //if(box.contains(i,j,k) || !boxed) {
      int nn = bod.nnodes(c);
                out << n+1 << " " << n+2 << " " << n+3 << boil::endl;
      if(nn>3)  out << n+1 << " " << n+3 << " " << n+4 << boil::endl;
      if(nn>4)  out << n+1 << " " << n+4 << " " << n+5 << boil::endl;
      if(nn>5)  out << n+1 << " " << n+5 << " " << n+6 << boil::endl;
      n+=nn;
  }
}

/******************************************************************************/
void PlotTEC::plot_tec_scalar(const Domain & dom, const Scalar & sca) {

  out << "# SCALAR VALUES" << boil::endl;

  int count=0;
  for(int k=boil::BW-sh; k<dom.nk()-boil::BW+sh; k++) 
    for(int j=boil::BW-sh; j<dom.nj()-boil::BW+sh; j++) 
      for(int i=boil::BW-sh; i<dom.ni()-boil::BW+sh; i++) {
        out << sca[i][j][k] << " ";       
        count++;
        if(count % 8 == 0) out << boil::endl;
      }
 
  /* end with new line (if needed) */
  if(count % 8 != 0) out << boil::endl;
}

/******************************************************************************/
void PlotTEC::plot_tec_scalarint(const Domain & dom, const ScalarInt & sca) {

  out << "# SCALAR VALUES" << boil::endl;

  int count=0;
  for(int k=boil::BW-sh; k<dom.nk()-boil::BW+sh; k++) 
    for(int j=boil::BW-sh; j<dom.nj()-boil::BW+sh; j++) 
      for(int i=boil::BW-sh; i<dom.ni()-boil::BW+sh; i++) {
        out << sca[i][j][k] << " ";       
        count++;
        if(count % 8 == 0) out << boil::endl;
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
    for(int k=boil::BW-sh; k<dom.nk()-boil::BW+sh; k++) 
      for(int j=boil::BW-sh; j<dom.nj()-boil::BW+sh; j++)
        for(int i=boil::BW-sh; i<dom.ni()-boil::BW+sh; i++) {
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
