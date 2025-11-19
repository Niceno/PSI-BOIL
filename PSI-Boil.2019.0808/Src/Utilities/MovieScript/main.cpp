#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include "../../Global/global_name_file.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  if(argc < 5 || argc > 9) {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] 
              << " <base_name> <first_step> <last_step> <stride>"  
              << " [fps=frames_per_second]"
              << " [aaf=anti-aliasing_factor]"
              << " [width=screen_width]"
              << " [name=movie_file_name]"
              << std::endl;
    return 0;
  }

  const char * bname  = argv[1];
  const char * first  = argv[2];       // first step           
  const char * last   = argv[3];       // last step           
  const int    stride = atoi(argv[4]); // stride

  /* set the default values for parameters */
  std::string width   ("720");
  std::string fps     ("12");
  std::string aaf     ("3");
  std::string name_avi("");
  std::string name_mcr("");

  /*===============================+
  |                                |
  |  analyze additional arguments  |
  |                                |
  +===============================*/
  if(argc > 5) {

    std::cout << "Processing additional arguments ..." << std::endl;

    /*-----------------------------------+
    |                                    |
    |  loop through aditional arguments  |
    |                                    |
    +-----------------------------------*/
    for(int arg=5; arg<argc; arg++) {
      std::string argument( argv[arg] ); 
      std::string option;
      std::string value;

      /*-----------------------------------+
      |  extract the option and the value  |
      +-----------------------------------*/
      for(int i=0; i < argument.length(); i++) {
        if( argument[i] == '=' ) {
          value  = argument.substr(i+1, argument.length()-1);
          option = argument.substr(0, i);
          break;
        } 
      }
//    std::cout << "option is: " << option << std::endl;
//    std::cout << "value  is: " << value  << std::endl;

      /*--------------------+
      |  frames per second  |
      +--------------------*/
      if( option == std::string("fps") ) { 
        fps = value;
      }

      /*--------+
      |  width  |
      +--------*/
      else if( option == std::string("width") ) { 
        width = value;
      }

      /*-------+
      |  name  |
      +-------*/
      else if( option == std::string("name") ) { 
        name_avi = value;
        name_mcr = value;
      }

      /*-----------------------+
      |  anti-aliasing factor  |
      +-----------------------*/
      else if( option == std::string("aaf") ) { 
        aaf = value;
      }

      /*---------------+
      |  unrecognized  |
      +---------------*/
      else {
        std::cout << "Unrecognized optoin " << option << "." << std::endl;
        exit(0);
      }
    } 
  }

  if(name_avi.length() == 0) {
    name_avi = std::string(bname)
             + std::string("_Step-") 
             + std::string(first) 
             + std::string("-") 
             + std::string(last) 
             + std::string("_Width-") 
             + std::string(width)
             + std::string("_FPS-") 
             + std::string(fps); 
  }

  /*---------------------------------------+ 
  |  write a little message on the screen  |
  +---------------------------------------*/
  std::cout << "Creating a Tecplot (TM) macro to create movie " << name_avi 
            << " with following parameters:" << std::endl;
  std::cout << "- frames per second   : " << fps   << std::endl;
  std::cout << "- width of the screen : " << width << std::endl;

  std::ifstream in;
  std::ofstream out;

  std::string line;
  std::string line_var;
  std::string var_list;

  /*-------------------------------+ 
  |  final settings for the names  |
  +-------------------------------*/
  name_mcr  = name_avi;  
  name_avi += std::string(".avi"); 
  name_mcr += std::string(".mcr"); 

  /*-------------------+
  |  << output header  |
  +-------------------*/
  out.open(name_mcr.c_str());
  out << "#!MC 1200" << std::endl;
  out << "# Created by Tecplot 360 build 12.2.0.9077" << std::endl;
  out << std::endl;

  /*=====================+
  |                      |
  |  loop through steps  |
  |                      |
  +=====================*/
  for(int step = atoi(first); step <= atoi(last); step += stride) {

    std::string current = name_file(bname, ".plt", step);

    /*------------------------------------------+
    |  << output commands for the current step  |
    +------------------------------------------*/
    out << "#------------------------"             << std::endl;
    out << "# Data for current frame "             << std::endl;
    out << "#------------------------"             << std::endl;
    out << "$!READDATASET \"" << current << "\""   << std::endl; 
    out << "  READDATAOPTION = NEW"                << std::endl;
    out << "  RESETSTYLE = NO"                     << std::endl;
    out << "  INCLUDETEXT = NO"                    << std::endl;
    out << "  INCLUDEGEOM = NO"                    << std::endl;
    out << "  INCLUDECUSTOMLABELS = NO"            << std::endl;
    out << "  VARLOADMODE = BYNAME"                << std::endl;
    out << "  ASSIGNSTRANDIDS = YES"               << std::endl;
    out << "  INITIALPLOTTYPE = CARTESIAN3D"       << std::endl;

    /*-------------+
    |              |
    |  first step  |
    |              |
    +-------------*/
    if(step == atoi(first)) {

      std::string current = name_file(bname, ".dat", step);

      /*--------------------+ 
      |  open current file  |
      +--------------------*/
      in.open(current.c_str());

      /* stop if file is not present */
      if( in.rdstate() != 0 ) {
        std::cout << "failed to open " << current << std::endl;
        std::cout << "exiting!" << std::endl;
        exit(0);
      }

      std::cout << "reading: " << current << std::endl;

      /*-----------------------+
      |  >> get variable info  | 
      +-----------------------*/

      /* skip first five lines */
      for(int l=0; l<5; l++) {
        getline(in, line);
//      std::cout << line << std::endl;
      }

      /* get line with variable info */
      getline(in, line_var);     /* " VARIABLES=... */
     
      /* extract variable names */
      int is, ie;
      for(int i=0; i < line_var.length(); i++) {
        if( line_var.at(i) == '=' ) {
          var_list = line_var.substr(i+1, line_var.length()-1);
          break;
        }
      }
      std::cout << "handling variables: " << var_list << std::endl;

      /*------------------------------------------+
      |  << output commands for the current step  |
      +------------------------------------------*/
      out << "  VARNAMELIST = '" << var_list << "'"  << std::endl;
      out << std::endl;

      /*--------------------------------------+
      |  << output settings fot the avi file  |
      +--------------------------------------*/
      out << "#-----------------------"                         << std::endl;
      out << "# Settings for AVI file "                         << std::endl;
      out << "#-----------------------"                         << std::endl;
      out << "$!EXPORTSETUP EXPORTFORMAT = AVI"                 << std::endl;
      out << "$!EXPORTSETUP IMAGEWIDTH = " << width             << std::endl;
      out << "$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES"   << std::endl;
      out << "$!EXPORTSETUP SUPERSAMPLEFACTOR = " << aaf        << std::endl;
      out << "$!EXPORTSETUP ANIMATIONSPEED = " << fps           << std::endl;
      out << "$!EXPORTSETUP AVICOMPRESSION = COLORPRESERVING"   << std::endl;
      out << "$!EXPORTSETUP EXPORTFNAME = '" << name_avi << "'" << std::endl;
      out << "$!EXPORTSTART"                                    << std::endl;
      out << "  EXPORTREGION = CURRENTFRAME"                    << std::endl;
      out << std::endl;

    /*------------------+
    |                   |
    |  all other steps  |
    |                   |
    +------------------*/
    } else {

      /*------------------------------------------+
      |  << output commands for the current step  |
      +------------------------------------------*/
      out << "  VARNAMELIST = '" << var_list << "'"  << std::endl;
      out << "$!EXPORTNEXTFRAME"                     << std::endl;
      out << std::endl;
    }

  }

  /*--------------------+
  |  << end the output  |
  +--------------------*/
  out << "#------------------------------" << std::endl;
  out << "# This line finalizes AVI file " << std::endl;
  out << "#------------------------------" << std::endl;
  out << "$!EXPORTFINISH "                 << std::endl;
}
