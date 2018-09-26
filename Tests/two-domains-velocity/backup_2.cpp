
    if((time.current_step()) % (nint)==0 ) {
      uvw_2  .save("uvw_2",   time.current_step());
      press_2.save("press_2", time.current_step());
      conc_2 .save("conc_2",  time.current_step());
      t_2  .save("t_2",       time.current_step());
      //nucl_2 .save("nucl",   time.current_step());
      //dmicro_2.save("dmicro",time.current_step());
    }
    if( boil::timer.current_min() > (wmin-30.0)
      || time.current_step()==time.total_steps()) {
      uvw_2  .save("uvw_2",   time.current_step());
      press_2.save("press_2", time.current_step());
      conc_2 .save("conc_2",  time.current_step());
      t_2    .save("t_2",     time.current_step());
      //nucl_2 .save("nucl_2",   time.current_step());
      //dmicro_2.save("dmicro_2",time.current_step());

      uvw_2  .rm("uvw_2",   ts);
      press_2.rm("press_2", ts);
      conc_2 .rm("conc_2",  ts);
      t_2    .rm("t_2",     ts);
      //nucl_2 .rm("nucl", ts);
      //dmicro_2.rm("dmicro", ts);
    }

