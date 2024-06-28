#include "pidcontrol.h"

/***************************************************************************//**
*  constructor of PIDcontrol class
*******************************************************************************/
real PIDcontrol::calculate(real setpoint,
                           real measured_value,
                           real current_time) {

  real dt = current_time - last_time_;
  real error = setpoint - measured_value;
  integral_ += error * dt;

  real derivative = 0.0;
  if (dt>0) {
    derivative = (error - prev_error_) / dt;
  }

  real output = kp_ * error + ki_ * integral_ + kd_ * derivative;

#if 0
  //std::cout<<"pid.cal "<<boil::cart.iam()<<" "<<last_time_<<" "<<dt<<"\n";
  //boil::oout<<"####pid.cal "<<boil::cart.iam()<<" "<<setpoint<<" "<<measured_value<<" "<<prev_error_<<"\n";
  //boil::oout<<"pid.cal "<<boil::cart.iam()<<" "<<error<<" "<<integral_<<" "<<derivative<<"\n";
#endif

  prev_error_ = error;
  last_time_ = current_time;

  return output;
}

