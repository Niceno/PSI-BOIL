#include "plot_tecmpi.h"
//#define READ_VEL_AVE  // velocity defined at face is average of cell center values
                        // If READ_VEL_AVE is not defined, then the cell center values
			// will be stored at face of uvw

std::string name_zone(int32_t zoneOwner){
  std::ostringstream oss;
  oss << "Proc_" << std::setw(5) << std::setfill('0') <<zoneOwner;
  return oss.str();
}

/******************************************************************************/
void PlotTECMPI::plot(Domain & dm, // couldn't make const out of "*this"
                   const char * nam, 
                   const int i, Times * t) {

  boil::timer.start("plotting");

  dom = & dm; // take it as a constant
  /* set domain size */
  plot_tecmpi_set_domain(dom);
  /* output file name */
  string fname = plot_tecmpi_fname(nam, i);
  boil::oout << "# Plotting: " << fname << boil::endl;
  /* variable names */
  std::vector<std::string> vnames;
  vnames.push_back("X");
  vnames.push_back("Y");
  vnames.push_back("Z");
  vnames.push_back("Rank");
  string variables = plot_tecmpi_variables(vnames);
  /* variable share */
  std::vector<int32_t> shareVarFromZone(numVars, 0); // No variable sharing for first zone output
  /* variable type */
  std::vector<int32_t> varTypes(numVars, FieldDataType_Float);  // set all variables float
  varTypes[3] = FieldDataType_Int32;  // variable no.4 is Rank = integer
  /* variable location */
  std::vector<int32_t> valueLocations(numVars, 1); // 0: cell center, 1: node point

  /* Open and initialize the file */
  plot_tecmpi_tecOpenInit(fname,variables);

  // Loop for zone
  for (int32_t zn = 0; zn < numZones; ++zn) {

    int32_t zone;
    int32_t zoneOwner = zn;  // owner (rank) of this zone

    // if zone-owner or commRank == 0, then crate file
    if (commRank == zn || commRank == mainRank) {
      /* Create zone and give map info */
      plot_tecmpi_tecCreateMap(varTypes,shareVarFromZone,valueLocations,zone,zoneOwner,i,t);
    }

    // if zone owner, then write data
    if (commRank == zn) {
      /* write data */
      int IVARNUM = 0;
      plot_tecmpi_write(dom,zone,IVARNUM);

      // write rank data
      std::vector<int32_t> rank(XDIM * YDIM * ZDIM, commRank);
      IVARNUM++;
      res = tecZoneVarWriteInt32Values(fileHandle, zone, IVARNUM, 0, rank.size(), &rank[0]);
    }
  }

  /* close file */
  res = tecFileWriterClose(&fileHandle);

  boil::timer.stop("plotting");
}

/******************************************************************************/
void PlotTECMPI::read(const char * nam,
                      const int i,
                      Times * t,
                      Vector * vec,
                      Scalar * sca,
                      Scalar * scb,
                      Scalar * scc,
                      Scalar * scd,
                      Scalar * sce,
                      Scalar * scf,
                      Scalar * scg,
                      Scalar * sch,
                      Scalar * sci) {
  dom = vec->domain(); // take it as a constant
  /* set domain size */
  plot_tecmpi_set_domain(dom);
  /* output file name */
  string fname = plot_tecmpi_fname(nam, i);
  boil::oout << "# Reading: " << fname << boil::endl;

  try {
    /* open file */
    res = tecFileReaderOpen(fname.c_str(), &fileHandle);

    /* variable names */
    int32_t numVars;
    res = tecDataSetGetNumVars(fileHandle, &numVars);
    boil::oout<<"# read:numVars= "<<numVars<<"\n";
    std::vector<std::string> vnames(numVars + 1);

    for (int32_t var = 1; var <= numVars; ++var) {
      char* name = NULL;
      res = tecVarGetName(fileHandle, var, &name);
      vnames[var] = name;
      //boil::oout<<" vnames,i= "<<vnames[var]<<" "<<var<<"\n";
    }

    // numZones
    res = tecDataSetGetNumZones(fileHandle, &numZones);
    boil::oout<<"# read:numZones= "<<numZones<<"\n";

    // solution time
    int32_t inputZone = 1;
    //real solutionTime;
    double solutionTime;
    res = tecZoneGetSolutionTime(fileHandle, inputZone, &solutionTime);
    t->current_time(solutionTime);
    boil::oout<<"# read:solutionTime= "<<real(solutionTime)<<"\n";

#ifdef READ_VEL_AVE
    /* reset vel */
    for_m(m)
      for_avmijk((*vec),m,i,j,k) {
        (*vec)[m][i][j][k]=0.0;
      }
#endif

    inputZone = commRank +1;
    for (int32_t var = 4; var <= numVars; ++var) {
      const int LEN = XDIM_C * YDIM_C * ZDIM_C;
      int64_t numValuesRead = 0;
      int64_t numValuesToRead = LEN;
      std::unique_ptr<float[]> values(new float[LEN]);  // receive

      /* read data */
      res = tecZoneVarGetFloatValues(fileHandle, inputZone, var, numValuesRead + 1,
                                     numValuesToRead, &values[0]);

#ifdef READ_VEL_AVE
      int icount[XDIM][YDIM][ZDIM];
      if(var==4) {
        Comp m = Comp::u();
        std::memset(icount, 0, sizeof(icount));
        for (int i = 0; i < XDIM_C; ++i)
          for (int j = 0; j < YDIM_C; ++j)
            for (int k = 0; k < ZDIM_C; ++k) {
              int index = (k * YDIM_C + j) * XDIM_C + i;
              (*vec)[m][i+BW  ][j+BW][k+BW] = values[index];
              (*vec)[m][i+BW+1][j+BW][k+BW] = values[index];
              icount[i  ][j][k]++;
              icount[i+1][j][k]++;
            }
        for (int i = 0; i < XDIM; ++i)
          for (int j = 0; j < YDIM; ++j)
            for (int k = 0; k < ZDIM; ++k) {
              if (icount[i][j][k] != 0)
                (*vec)[m][i+BW][j+BW][k+BW] /= icount[i][j][k];
            }
      }
      if(var==5) {
        Comp m = Comp::v();
        std::memset(icount, 0, sizeof(icount));
        for (int i = 0; i < XDIM_C; ++i)
          for (int j = 0; j < YDIM_C; ++j)
            for (int k = 0; k < ZDIM_C; ++k) {
              int index = (k * YDIM_C + j) * XDIM_C + i;
              (*vec)[m][i+BW][j+BW  ][k+BW] = values[index];
              (*vec)[m][i+BW][j+BW+1][k+BW] = values[index];
              icount[i][j  ][k]++;
              icount[i][j+1][k]++;
            }
        for (int i = 0; i < XDIM; ++i)
          for (int j = 0; j < YDIM; ++j)
            for (int k = 0; k < ZDIM; ++k) {
              if (icount[i][j][k] != 0)
                (*vec)[m][i+BW][j+BW][k+BW] /= icount[i][j][k];
            }
      }
      if(var==6) {
        Comp m = Comp::w();
        std::memset(icount, 0, sizeof(icount));
        for (int i = 0; i < XDIM_C; ++i)
          for (int j = 0; j < YDIM_C; ++j)
            for (int k = 0; k < ZDIM_C; ++k) {
              int index = (k * YDIM_C + j) * XDIM_C + i;
              (*vec)[m][i+BW][j+BW][k+BW  ] = values[index];
              (*vec)[m][i+BW][j+BW][k+BW+1] = values[index];
              icount[i][j][k  ]++;
              icount[i][j][k+1]++;
            }
        for (int i = 0; i < XDIM; ++i)
          for (int j = 0; j < YDIM; ++j)
            for (int k = 0; k < ZDIM; ++k) {
              if (icount[i][j][k] != 0)
                (*vec)[m][i+BW][j+BW][k+BW] /= icount[i][j][k];
            }
      }
#endif

#ifndef READ_VEL_AVE
      if(var==4 ) copy_valCell(values,vec,Comp::u()); 
      if(var==5 ) copy_valCell(values,vec,Comp::v()); 
      if(var==6 ) copy_valCell(values,vec,Comp::w()); 
#endif
      if(var==7 ) copy_valCell(values,sca);
      if(var==8 ) copy_valCell(values,scb);
      if(var==9 ) copy_valCell(values,scc);
      if(var==10) copy_valCell(values,scd);
      if(var==11) copy_valCell(values,sce);
      if(var==12) copy_valCell(values,scf);
      if(var==13) copy_valCell(values,scg);
      if(var==14) copy_valCell(values,scg);
      if(var==15) copy_valCell(values,scg);
    }

    /* close file */
    res = tecFileReaderClose(&fileHandle);
  } catch (std::runtime_error const& e) {
    std::cerr << "Error: " << e.what() << "Proc=" <<commRank<< std::endl;
    exit(0);
  }

}
/******************************************************************************/
void PlotTECMPI::read(const char * nam,
                      const int i,
                      Times * t,
                      Scalar * sca,
                      Scalar * scb,
                      Scalar * scc,
                      Scalar * scd,
                      Scalar * sce,
                      Scalar * scf,
                      Scalar * scg,
                      Scalar * sch,
                      Scalar * sci) {
  dom = sca->domain(); // take it as a constant
  /* set domain size */
  plot_tecmpi_set_domain(dom);
  /* output file name */
  string fname = plot_tecmpi_fname(nam, i);
  boil::oout << "# Reading: " << fname << boil::endl;

  try {
    /* open file */
    res = tecFileReaderOpen(fname.c_str(), &fileHandle);

    /* variable names */
    int32_t numVars;
    res = tecDataSetGetNumVars(fileHandle, &numVars);
    boil::oout<<"# read:numVars= "<<numVars<<"\n";
    std::vector<std::string> vnames(numVars + 1);

    for (int32_t var = 1; var <= numVars; ++var) {
      char* name = NULL;
      res = tecVarGetName(fileHandle, var, &name);
      vnames[var] = name;
      //boil::oout<<" vnames,i= "<<vnames[var]<<" "<<var<<"\n";
    }

    // numZones
    res = tecDataSetGetNumZones(fileHandle, &numZones);
    boil::oout<<"# read:numZones= "<<numZones<<"\n";

    // solution time
    int32_t inputZone = 1;
    //real solutionTime;
    double solutionTime;
    res = tecZoneGetSolutionTime(fileHandle, inputZone, &solutionTime);
    t->current_time(solutionTime);
    boil::oout<<"# read:solutionTime= "<<real(solutionTime)<<"\n";

    inputZone = commRank +1;
    for (int32_t var = 4; var <= numVars; ++var) {
      const int LEN = XDIM_C * YDIM_C * ZDIM_C;
      int64_t numValuesRead = 0;
      int64_t numValuesToRead = LEN;
      std::unique_ptr<float[]> values(new float[LEN]);  // receive

      /* read data */
      res = tecZoneVarGetFloatValues(fileHandle, inputZone, var, numValuesRead + 1,
                                     numValuesToRead, &values[0]);

      if(var==4 ) copy_valCell(values,sca);
      if(var==5 ) copy_valCell(values,scb);
      if(var==6 ) copy_valCell(values,scc);
      if(var==7) copy_valCell(values,scd);
      if(var==8) copy_valCell(values,sce);
      if(var==9) copy_valCell(values,scf);
      if(var==10) copy_valCell(values,scg);
      if(var==11) copy_valCell(values,scg);
      if(var==12) copy_valCell(values,scg);
    }

    /* close file */
    res = tecFileReaderClose(&fileHandle);
  } catch (std::runtime_error const& e) {
    std::cerr << "Error: " << e.what() << "Proc=" <<commRank<< std::endl;
    exit(0);
  }
}
/******************************************************************************/
void PlotTECMPI::plot(const char * nam,
                      const int i,
                      Times * t,
                      const Scalar * sca,
                      const Scalar * scb,
                      const Scalar * scc,
                      const Scalar * scd,
                      const Scalar * sce,
                      const Scalar * scf,
                      const Scalar * scg,
                      const Scalar * sch,
                      const Scalar * sci) {
  boil::timer.start("plotting");

  dom = sca->domain(); // take it as a constant
  /* set domain size */
  plot_tecmpi_set_domain(dom);
  /* output file name */
  string fname = plot_tecmpi_fname(nam, i);
  boil::oout << "# Plotting: " << fname << boil::endl;
  /* variable names */
  std::vector<std::string> vnames;
  vnames.push_back("X");
  vnames.push_back("Y");
  vnames.push_back("Z");
  vnames.push_back("A"); if(sca->name().length() > 0) vnames [3] = sca->name();
  if(scb!=NULL) {vnames.push_back("B"); if(scb->name().length() > 0) vnames [4] = scb->name();}
  if(scc!=NULL) {vnames.push_back("C"); if(scc->name().length() > 0) vnames [5] = scc->name();}
  if(scd!=NULL) {vnames.push_back("D"); if(scd->name().length() > 0) vnames [6] = scd->name();}
  if(sce!=NULL) {vnames.push_back("E"); if(sce->name().length() > 0) vnames [7] = sce->name();}
  if(scf!=NULL) {vnames.push_back("F"); if(scf->name().length() > 0) vnames [8] = scf->name();}
  if(scg!=NULL) {vnames.push_back("G"); if(scg->name().length() > 0) vnames [9] = scg->name();}
  if(sch!=NULL) {vnames.push_back("H"); if(sch->name().length() > 0) vnames [10] = sch->name();}
  if(sci!=NULL) {vnames.push_back("I"); if(sci->name().length() > 0) vnames [11] = sci->name();}
  string variables = plot_tecmpi_variables(vnames);
  /* variable share */
  std::vector<int32_t> shareVarFromZone(numVars, 0); // No variable sharing for first zone output
  /* variable type */
  std::vector<int32_t> varTypes(numVars, FieldDataType_Float);  // set all variables float
  /* value location */
  std::vector<int32_t> valueLocations(numVars, 0); // 0: cell center, 1: node point
  valueLocations[0] = 1;  // X node point
  valueLocations[1] = 1;  // Y node point
  valueLocations[2] = 1;  // Z node point

  /* Open and initialize the file */
  plot_tecmpi_tecOpenInit(fname,variables);

  // Loop for zone
  for (int32_t zn = 0; zn < numZones; ++zn) {

    int32_t zone;
    int32_t zoneOwner = zn;  // owner (rank) of this zone

    // if zone-owner or commRank == 0, then crate file
    if (commRank == zn || commRank == mainRank) {
      /* Create zone and give map info */
      plot_tecmpi_tecCreateMap(varTypes,shareVarFromZone,valueLocations,zone,zoneOwner,i,t);
    }

    // if zone owner, then write data
    if (commRank == zn) {
      int IVARNUM = 0;
      plot_tecmpi_write(dom,zone,IVARNUM);
      plot_tecmpi_write(sca,zone,IVARNUM);
      if(scb!=NULL) plot_tecmpi_write(scb,zone,IVARNUM);
      if(scc!=NULL) plot_tecmpi_write(scc,zone,IVARNUM);
      if(scd!=NULL) plot_tecmpi_write(scd,zone,IVARNUM);
      if(sce!=NULL) plot_tecmpi_write(sce,zone,IVARNUM);
      if(scf!=NULL) plot_tecmpi_write(scf,zone,IVARNUM);
      if(scg!=NULL) plot_tecmpi_write(scg,zone,IVARNUM);
      if(sch!=NULL) plot_tecmpi_write(sch,zone,IVARNUM);
      if(sci!=NULL) plot_tecmpi_write(sci,zone,IVARNUM);
    }
  }

  /* close file */
  res = tecFileWriterClose(&fileHandle);

  boil::timer.stop("plotting");
}
/******************************************************************************/
void PlotTECMPI::plot(const char * nam,
                      const int i,
                      Times * t,
                      const Vector * vec,
                      const Scalar * sca,
                      const Scalar * scb,
                      const Scalar * scc,
                      const Scalar * scd,
                      const Scalar * sce,
                      const Scalar * scf,
                      const Scalar * scg,
                      const Scalar * sch,
                      const Scalar * sci) {
  boil::timer.start("plotting");

  dom = vec->domain(); // take it as a constant
  /* set domain size */
  plot_tecmpi_set_domain(dom);
  /* output file name */
  string fname = plot_tecmpi_fname(nam, i);
  boil::oout << "# Plotting: " << fname << boil::endl;
  /* variable names */
  std::vector<std::string> vnames;
  vnames.push_back("X");
  vnames.push_back("Y");
  vnames.push_back("Z");
  vnames.push_back("U");
  vnames.push_back("V");
  vnames.push_back("W");
  if(sca!=NULL) {vnames.push_back("A"); if(sca->name().length() > 0) vnames [6] = sca->name();}
  if(scb!=NULL) {vnames.push_back("B"); if(scb->name().length() > 0) vnames [7] = scb->name();}
  if(scc!=NULL) {vnames.push_back("C"); if(scc->name().length() > 0) vnames [8] = scc->name();}
  if(scd!=NULL) {vnames.push_back("D"); if(scd->name().length() > 0) vnames [9] = scd->name();}
  if(sce!=NULL) {vnames.push_back("E"); if(sce->name().length() > 0) vnames [10] = scd->name();}
  if(scf!=NULL) {vnames.push_back("F"); if(scf->name().length() > 0) vnames [11] = scd->name();}
  if(scg!=NULL) {vnames.push_back("G"); if(scg->name().length() > 0) vnames [12] = scd->name();}
  if(sch!=NULL) {vnames.push_back("H"); if(sch->name().length() > 0) vnames [13] = scd->name();}
  if(sci!=NULL) {vnames.push_back("I"); if(sci->name().length() > 0) vnames [14] = scd->name();}
  string variables = plot_tecmpi_variables(vnames);
  /* variable share */
  std::vector<int32_t> shareVarFromZone(numVars, 0); // No variable sharing for first zone output
  /* variable type */
  std::vector<int32_t> varTypes(numVars, FieldDataType_Float);  // set all variables float
  /* value location */
  std::vector<int32_t> valueLocations(numVars, 0); // 0: cell center, 1: node point
  valueLocations[0] = 1;  // X node point
  valueLocations[1] = 1;  // Y node point
  valueLocations[2] = 1;  // Z node point

  /* Open and initialize the file */
  plot_tecmpi_tecOpenInit(fname,variables);

  // Loop for zone
  for (int32_t zn = 0; zn < numZones; ++zn) {

    int32_t zone;
    int32_t zoneOwner = zn;  // owner (rank) of this zone

    // if zone-owner or commRank == 0, then crate file
    if (commRank == zn || commRank == mainRank) {
      /* Create zone and give map info */
      plot_tecmpi_tecCreateMap(varTypes,shareVarFromZone,valueLocations,zone,zoneOwner,i,t);
    }

    // if zone owner, then write data
    if (commRank == zn) {
      int IVARNUM = 0;
      plot_tecmpi_write(dom,zone,IVARNUM);
      plot_tecmpi_write(vec,zone,IVARNUM);
      if(sca!=NULL) plot_tecmpi_write(sca,zone,IVARNUM);
      if(scb!=NULL) plot_tecmpi_write(scb,zone,IVARNUM);
      if(scc!=NULL) plot_tecmpi_write(scc,zone,IVARNUM);
      if(scd!=NULL) plot_tecmpi_write(scd,zone,IVARNUM);
      if(sce!=NULL) plot_tecmpi_write(sce,zone,IVARNUM);
      if(scf!=NULL) plot_tecmpi_write(scf,zone,IVARNUM);
      if(scg!=NULL) plot_tecmpi_write(scg,zone,IVARNUM);
      if(sch!=NULL) plot_tecmpi_write(sch,zone,IVARNUM);
      if(sci!=NULL) plot_tecmpi_write(sci,zone,IVARNUM);
    }
  }

  /* close file */
  res = tecFileWriterClose(&fileHandle);

  boil::timer.stop("plotting");
}
/******************************************************************************/
void PlotTECMPI::plot(const Pathline & pl, 
                   const char * nam, 
                   const int i, Times * t) {
  if( boil::cart.iam()==0 ) {
    /* open the result file */
    std::string name = name_file(nam, ".dat", i);
    boil::oout << "# Plotting: "<<name<<"\n";
    out.open(name.c_str());

    // visit cannot read the next line.
    //out << "#Number_of_variables= "<<6+pl.nval()<<"\n";
    out << "VARIABLES= X Y Z U V W ID ";
    if(pl.dia_den()) {
      out << "Diameter Density ";
    }
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
          << pl.particles[ip].w()<<" "
          << pl.particles[ip].id()<<" ";
      if(pl.dia_den()) {
        out << pl.particles[ip].diameter()<<" "
            << pl.particles[ip].density()<<" ";
      }
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
std::string PlotTECMPI::plot_tecmpi_fname(const char * nam, const int i) {

  /* file name extension */
  std::string name = name_file(nam, ".szplt", i); // i: time step
  return name;
}

/******************************************************************************/
void PlotTECMPI::plot_tecmpi_set_domain(const Domain * dom){

  /* set array size for decomposed domain */
  XDIM = dom->ni() -2*BW +1;
  YDIM = dom->nj() -2*BW +1;
  ZDIM = dom->nk() -2*BW +1;

  XDIM_C = XDIM -1;
  YDIM_C = YDIM -1;
  ZDIM_C = ZDIM -1;
} 

/******************************************************************************/
std::string PlotTECMPI::plot_tecmpi_variables
                        (const std::vector<std::string>& vnames) {
  // set numVars
  numVars = vnames.size();
  // Use a string stream to concatenate the vector elements
  std::ostringstream oss;
  for (size_t i = 0; i < vnames.size(); ++i) {
    oss << vnames[i];
    if (i < vnames.size() - 1) {
        oss << " "; // Add a space between elements
    }
  }
  //boil::oout<<"plot_tecmpi_variables:variables= "<<oss.str()
  //          <<"\nnumVars= "<<numVars<<"\n";

  // Return the concatenated result as a string
  return oss.str();
}

/******************************************************************************/
void PlotTECMPI::plot_tecmpi_tecOpenInit(const string & fname,
                                         const string & variables){
  /* Open the file and write the tecplot datafile header information */
  res = tecFileWriterOpen(
        fname.c_str(),
        "IJK Ordered Zone",
        variables.c_str(),
        FILEFORMAT_SZL,
        FILETYPE_FULL,
        FieldDataType_Float,
        NULL,
        &fileHandle);

  /* Initialize tecMPI done by all processes */
  res = tecMPIInitialize(fileHandle, mpiComm, mainRank);
}

/******************************************************************************/
void PlotTECMPI::plot_tecmpi_tecCreateMap(vector<int32_t> & varTypes,
                                          vector<int32_t> & shareVarFromZone,
                                          vector<int32_t> & valueLocations,
                                          int32_t & zone, 
                                          int32_t & zoneOwner,
                                          const int & strandID,
                                          const Times * t) {

  res = tecZoneCreateIJK(
          fileHandle,
	  name_zone(zoneOwner).c_str(),  // Zone name 
          XDIM,
          YDIM,
          ZDIM,
          &varTypes[0],
          &shareVarFromZone[0],
          &valueLocations[0],
          NULL, // passiveVarList
          0,    // shareFaceNeighborsFromZone,
          0,    // numFaceConnections,
          0,    // faceNeighborMode,
          &zone);

  res = tecZoneMapPartitionsToMPIRanks(
          fileHandle,
          zone,
          1, // numPartitions,
          &zoneOwner);

  if(t != NULL){
    res = tecZoneSetUnsteadyOptions(
          fileHandle,
          zone,
          t->current_time(),
          strandID);
  }
}

/******************************************************************************/
void PlotTECMPI::plot_tecmpi_write(const Domain * dom, int32_t & zone,
                                   int & IVARNUM) {

  const int LEN = XDIM * YDIM * ZDIM;
  std::unique_ptr<float[]> xyz(new float[LEN]);

  // X
  //std::unique_ptr<float[]> xyz(new float[LEN]);
  // Calculate the zone variables
  for (int i = 0; i < XDIM; ++i)
    for (int j = 0; j < YDIM; ++j)
      for (int k = 0; k < ZDIM; ++k) { 
        int index = (k * YDIM + j) * XDIM + i;
        xyz[index] = dom->xn(i+BW);
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &xyz[0]);

  // Y
  //xyz.reset(new float[LEN]);
  // Calculate the zone variables
  for (int i = 0; i < XDIM; ++i)
    for (int j = 0; j < YDIM; ++j)
      for (int k = 0; k < ZDIM; ++k) { 
        int index = (k * YDIM + j) * XDIM + i;
        xyz[index] = dom->yn(j+BW);
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &xyz[0]);

  // Z
  //xyz.reset(new float[LEN]);
  // Calculate the zone variables
  for (int i = 0; i < XDIM; ++i)
    for (int j = 0; j < YDIM; ++j)
      for (int k = 0; k < ZDIM; ++k) {
        int index = (k * YDIM + j) * XDIM + i;
        xyz[index] = dom->zn(k+BW);
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &xyz[0]);
}

/******************************************************************************/
void PlotTECMPI::plot_tecmpi_write(const Vector * vec, int32_t & zone,
                                   int & IVARNUM) {

  const int LEN = XDIM_C * YDIM_C * ZDIM_C;
  std::unique_ptr<float[]> uvw(new float[LEN]);

  // U
  // Calculate the zone variables
  Comp m = Comp::u();
  for (int i = 0; i < XDIM_C; ++i)
    for (int j = 0; j < YDIM_C; ++j)
      for (int k = 0; k < ZDIM_C; ++k) {
        int index = (k * YDIM_C + j) * XDIM_C + i;
        uvw[index] = 0.5 * ((*vec)[m][i+BW][j+BW][k+BW] + (*vec)[m][i+1+BW][j+BW][k+BW]);
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &uvw[0]);

  // V
  // Calculate the zone variables
  m = Comp::v(); 
  for (int i = 0; i < XDIM_C; ++i)
    for (int j = 0; j < YDIM_C; ++j)
      for (int k = 0; k < ZDIM_C; ++k) {
        int index = (k * YDIM_C + j) * XDIM_C + i;
        uvw[index] = 0.5 * ((*vec)[m][i+BW][j+BW][k+BW] + (*vec)[m][i+BW][j+1+BW][k+BW]);
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &uvw[0]);

  // W
  // Calculate the zone variables
  m = Comp::w();
  for (int i = 0; i < XDIM_C; ++i)
    for (int j = 0; j < YDIM_C; ++j)
      for (int k = 0; k < ZDIM_C; ++k) {
        int index = (k * YDIM_C + j) * XDIM_C + i;
        uvw[index] = 0.5 * ((*vec)[m][i+BW][j+BW][k+BW] + (*vec)[m][i+BW][j+BW][k+1+BW]);
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &uvw[0]);
}

/******************************************************************************/
void PlotTECMPI::plot_tecmpi_write(const Scalar * s, int32_t & zone,
                                   int & IVARNUM) {

  const int BW = boil::BW;
  const int LEN = XDIM_C * YDIM_C * ZDIM_C;
  std::unique_ptr<float[]> val(new float[LEN]);

  // Calculate the zone variables
  for (int i = 0; i < XDIM_C; ++i)
    for (int j = 0; j < YDIM_C; ++j)
      for (int k = 0; k < ZDIM_C; ++k) {
        int index = (k * YDIM_C + j) * XDIM_C + i;
        val[index] = (*s)[i+BW][j+BW][k+BW];
      }
  IVARNUM++;
  //std::cout<<"plot_tecmpi_write:Rank= "<<commRank<<" IVARNUM "<<IVARNUM<<" zone "<<zone<<"\n";
  res = tecZoneVarWriteFloatValues(fileHandle, zone, IVARNUM, 0, LEN, &val[0]);

}

/******************************************************************************/
void PlotTECMPI::copy_valCell(const std::unique_ptr<float[]>& values, Scalar *s) {
  for (int i = 0; i < XDIM_C; ++i)
    for (int j = 0; j < YDIM_C; ++j)
      for (int k = 0; k < ZDIM_C; ++k) {
        int index = (k * YDIM_C + j) * XDIM_C + i;
        (*s)[i+BW][j+BW][k+BW] = values[index];
      }
}

/******************************************************************************/
void PlotTECMPI::copy_valCell(const std::unique_ptr<float[]>& values,
                              Vector *v, const Comp &m) {
  for (int i = 0; i < XDIM_C; ++i)
    for (int j = 0; j < YDIM_C; ++j)
      for (int k = 0; k < ZDIM_C; ++k) {
        int index = (k * YDIM_C + j) * XDIM_C + i;
        (*v)[m][i+BW][j+BW][k+BW] = values[index];
      }
  boil::oout<<"copy_valCell: "<<values[0]<<" "<<(*v)[m][BW][BW][BW]<<"\n";
}
