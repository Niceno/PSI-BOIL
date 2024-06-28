#include "pidcontrol.h"

/***************************************************************************//**
*  constructor of PIDcontrol class
*
*  kp: coefficient for proportional term
*  ki: coefficient for integral term
*  kd: coefficient for delivertive term
*******************************************************************************/
PIDcontrol::PIDcontrol(real kp, real ki, real kd)
        : kp_(kp), ki_(ki), kd_(kd),
          prev_error_(0.0), integral_(0.0), last_time_(0.0) 
{
}

