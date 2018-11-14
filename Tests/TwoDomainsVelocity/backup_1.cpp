
    if((time.current_step()) % (nint)==0 ) {
      uvw_1  .save("uvw_1",   time.current_step());
      press_1.save("press_1", time.current_step());
      conc_1 .save("conc_1",  time.current_step());
      t_1  .save("t_1",       time.current_step());
      //nucl_1 .save("nucl",   time.current_step());
      //dmicro_1.save("dmicro",time.current_step());
    }
    if( boil::timer.current_min() > (wmin-30.0)
      || time.current_step()==time.total_steps()) {
      uvw_1  .save("uvw_1",   time.current_step());
      press_1.save("press_1", time.current_step());
      conc_1 .save("conc_1",  time.current_step());
      t_1    .save("t_1",     time.current_step());
      //nucl_1 .save("nucl_1",   time.current_step());
      //dmicro_1.save("dmicro_1",time.current_step());

      uvw_1  .rm("uvw_1",   ts);
      press_1.rm("press_1", ts);
      conc_1 .rm("conc_1",  ts);
      t_1    .rm("t_1",     ts);
      //nucl_1 .rm("nucl", ts);
      //dmicro_1.rm("dmicro", ts);
    }

