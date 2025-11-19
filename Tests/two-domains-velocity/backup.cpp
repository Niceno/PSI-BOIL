
    if((time.current_step()) % (nint)==0 ) {
      if( boil::cart.iam()==0) {
        std::fstream output;
        std::stringstream ss;
        ss <<"time-"<<time.current_step()<<".txt";
        std::string fname = ss.str();
        int len = fname.length();
        char * cfname = new char[len+1];
        memcpy(cfname, fname.c_str(), len+1);
        output << std::setprecision(16);
        output.open(cfname, std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output.close();
      }
    }
    if( boil::timer.current_min() > (wmin-30.0)
      || time.current_step()==time.total_steps()) {
      std::fstream output;
      output << std::setprecision(16);
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      output.open("run.txt", std::ios::out);
      output << 0 << boil::endl;
      output.close();

      boil::timer.stop();
      boil::timer.report();
      exit(0);

    }

