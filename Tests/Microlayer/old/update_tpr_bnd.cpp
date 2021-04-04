    /* reset boundary condition for temperature */
    if(NZsol==0&&(case_flag==2||case_flag==4)) {
      std::vector<real> C0 = {12.56791608260579,
                              -0.38231422424431977,
                              -2.5843730362760384,
                              1.0982624468028859,
                              -0.1961483629504791};
      std::vector<real> C1 = {8.302991250503421,
                              13.176340482203388,
                              -15.81871438800867,
                              5.836534846472568,
                              -0.6514010904922061};

      std::vector<real> Cinter = C0;
      for(int i(0); i<Cinter.size();++i) {
        Cinter[i] += time.current_time()/0.21e-3 * (C1[i]-C0[i]);
      }

      std::vector<std::ostringstream> sci(Cinter.size());
      std::vector<std::string> scistr(Cinter.size());
      std::ostringstream fullstr;
      boil::oout<<"tpr_bnd_update= "<<time.current_time()<<" ";
      for(int i(0); i<Cinter.size();++i) {
        sci[i]<<Cinter[i];
        scistr[i] = sci[i].str();
        boil::oout<<scistr[i]<<" ";
        fullstr<<"("<<scistr[i]<<")*(((x/1e-3)^2+(y/1e-3)^2)^0.5)^"<<i<<"+";
      }
      boil::oout<<boil::endl;
      fullstr<<"0.";
      //boil::oout<<"eq= "<<fullstr.str()<<boil::endl;
      char *eqtpr = new char[fullstr.str().length()+1];
      std::strcpy(eqtpr, fullstr.str().c_str());

      for(auto l : tpr.levels) {
       l->bc().modify( BndCnd( Dir::kmin(), BndType::dirichlet(), eqtpr) );
       l->bnd_update();
     }
    }
