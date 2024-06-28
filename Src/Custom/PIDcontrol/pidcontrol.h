#ifndef PIDCONTROL_H
#define PIDCONTROL_H
#include <iostream>
#include <chrono>
#include <thread>
#include "pidcontrol.h"
#include "../../Global/global_precision.h"
#include "../../Global/global_name_file.h"
#include "../../Parallel/communicator.h"
#include "../../Parallel/mpi_macros.h"
#include "../../Parallel/Out/print.h"

///////////////////////////////////////////////////
//                                               //
//  PID controler                                //
//  Proportional–integral–derivative controller  //
//                                               //
///////////////////////////////////////////////////
class PIDcontrol {
public:
    PIDcontrol(real kp, real ki, real kd);

    real calculate(real setp, real meas, real time);
    void load(const char * nm, const int it);
    void load(std::ifstream & in);
    void rm(const char * nm, const int it);
    void save(const char * nm, const int it);
    void save(std::ofstream & out);


private:
    real kp_;
    real ki_;
    real kd_;
    real prev_error_;
    real integral_;
    real last_time_;
};
#endif

