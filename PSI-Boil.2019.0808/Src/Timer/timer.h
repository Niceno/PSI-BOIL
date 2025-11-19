/***************************************************************************//**
*  \brief Global class for timing different parts of the code.
*
*  Designed for profiling programs. Essentially, contains two subroutines:
*  "start" and "stop". If they are are invoked without arguments, they will
*  time the execution of the program. If invoked with a name, they will time
*  only part of the code (tipically a function, but parts of the loops may
*  also be timed). 
*
*  In addition to these two essential functions ("start" and "stop"), it 
*  also contains routines to report the meassured times ("report").
*  Typically, calling "report" makes sense only at the end of the program.
*
*  \warning
*  All local timings, started with "start(name)", must be stopped with
*  a corresponding "stop" name. Unfortunatelly, one has to take special care
*  to fulfill this, as class does not check it.
*******************************************************************************/

#ifndef TIMER_H
#define TIMER_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include <ctime>
#include <vector>

#include "../Global/global_precision.h"
#include "../Parallel/Out/out.h"
#include "../Parallel/Out/print.h"
#include "../Global/global_minmax.h"

/////////////
//         //
//  Timer  //
//         //
/////////////
class Timer {
  public:
    Timer();
    void start(); /* for global */
    void stop();  /* fot global */
    void start(const std::string fname); /* for local */
    void stop(const std::string fname);  /* for local */
                                      
    void report();
    void report_separate();

    real current_sec();
    real current_min();
    real current_hour();

  private:
    real elapsed_s();
    bool               is_running;
    int                i_stopped;
    std::clock_t       start_time;
    std::clock_t       total_time;
    std::string        name;
    std::vector<Timer> sub; /* local (subroutine) timers */
};

/*---------------+
|  global timer  |
+---------------*/
namespace boil {
  extern Timer timer;
}

#endif
