#include "custom.h"

namespace boil {
  void test_irun(const std::string & testfile) {
    std::fstream input;
    int irun(0);

    /* get irun value */
    if(boil::cart.iam()==0) {
      input.open(testfile, std::ios::in);
      if( !input.fail() ) {
        input >> irun;
        boil::oout<<"read irun. irun= "<<irun<<boil::endl;
      }
      input.close();
    }
    boil::cart.sum_int(&irun);

    /* exit? */
    if(irun==1) {
      boil::oout<<"exit job due to irun=1"<<boil::endl;
      exit(0);
    }

    return;
  }

  void set_irun(const int val, const std::string & testfile) {
    if(boil::cart.iam()==0) {
      std::fstream output;
      output.open(testfile, std::ios::out);
      output << val << boil::endl;
      output.close();
    }
  }

} /* namespace */
