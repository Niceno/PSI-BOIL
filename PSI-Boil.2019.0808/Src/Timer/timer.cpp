#include "timer.h" 

/*---------------+
|  global timer  |
+---------------*/
namespace boil {
  Timer timer;
}

/***************************************************************************//**
*  Basic constructor. Does not start the timing at the time it is created. 
*******************************************************************************/
Timer :: Timer() {
  is_running = false;
  i_stopped  = -1;
  total_time = 0;
  name       = "";
}

/***************************************************************************//**
*  This starts the global (unnamed) counter.   
*******************************************************************************/
void Timer :: start() {
  if(is_running)
    boil::oout << "timer " << name << " is allready running !" << boil::endl;
  else {
    is_running = true; 
    start_time = clock();
  }
}

/***************************************************************************//**
*  Stops the global (unnamed) counter.      
*******************************************************************************/
void Timer :: stop() {
  if(!is_running)
    boil::oout << "timer " << name << " is not running !" << boil::endl;
  else {
    total_time += (clock() - start_time);
    is_running = false;
  }
}

/***************************************************************************//**
*  Creates and starts local (named) timer if it has not been created already.
*  If the local timer has been created, then it only starts it. 
*  In either case, stops other local timer.
********************************************************************************/
void Timer :: start(const std::string fname) {

  /* find it */
  int it=-1;
  for(int i=0; i<sub.size(); i++) {
    if( sub[i].name == fname ) {
      it=i;
      break;
    }
  }

  /* if it doesn't exist, create it */
  if(it == -1) {
    Timer newt;          /* new timer */
    newt.name = fname;   /* set its name */
    sub.push_back(newt);
    it = sub.size()-1;
  }

  /* stop the others, storing the stopped one */
  sub[it].i_stopped = -1;
  for(int i=0; i<sub.size(); i++) {
    if( sub[i].is_running ) {
      sub[i].stop();
      sub[it].i_stopped = i;
      break;
    }
  }

  /* start it */
  sub[it].start();
}

/***************************************************************************//**
*  Stops the local (subroutine) timer, and restarts the one you stopped. 
*******************************************************************************/
void Timer :: stop(const std::string fname) {

  /* if present stop it ... */
  for(int it=0; it<sub.size(); it++) {
    if( sub[it].name == fname ) { /* found it */

      sub[it].stop(); /* stop it */

      /* ... and start the one it stopped */
      if(sub[it].i_stopped != -1) {
        sub[sub[it].i_stopped].start();
        sub[it].i_stopped = -1;
      }
      return;
    }
  }

  /* if not present, give a warning message */
  boil::oout << "timer " << fname << " isn't running !" << boil::endl;
}

/***************************************************************************//**
*  Returns the current time in seconds.
*******************************************************************************/
real Timer :: current_sec() {

  real cur=0.0;

  if(!is_running)
    boil::oout << "timer " << name << " is not running !" << boil::endl;
  else 
    cur = (clock() - start_time) / (real)CLOCKS_PER_SEC;

  boil::cart.average_real(& cur);

  return cur;
}

/***************************************************************************//**
*  Returns the current time in minutes.
*******************************************************************************/
real Timer :: current_min() {

  return current_sec() / 60.0;
}

/***************************************************************************//**
*  Returns the current time in minutes.
*******************************************************************************/
real Timer :: current_hour() {

  return current_min() / 60.0;
}

/******************************************************************************/
real Timer :: elapsed_s() {
  real e = (real)total_time / (real)CLOCKS_PER_SEC;
  return e;
}

/***************************************************************************//**
*  Reports the time spent in the program. 
*******************************************************************************/
void Timer :: report_separate() {

  /* total time */
  real total = elapsed_s();
  boil::aout << "+==========================================" << boil::endl;
  boil::aout << "| Total execution time: " << total << " [s]" << boil::endl;
  boil::aout << "+------------------------------------------" << boil::endl;

  if(total < 0.01) return;

  std::string elsew = "elsewhere";

  int longest(0);
  for(int i=0; i<sub.size(); i++) 
    longest=boil::maxi(longest, sub[i].name.length());
  for(int i=0; i<sub.size(); i++) 
    sub[i].name.resize(longest, ' ');
  elsew.resize(longest+3, ' '); /* 3 is length of "in " */

  real acum = 0.0; 
  for(int i=0; i<sub.size(); i++) {
    real here = sub[i].elapsed_s();
    boil::aout << "| Time spent in " << sub[i].name << ": " 
               << here << " [s]    (" 
               << here * 100.0 / total << "%)" << boil::endl;
    acum += here;
  }
  if(sub.size()) {
    boil::aout << "| Time spent " << elsew << ": " 
               << total - acum << " [s]    ("
               << (total - acum) * 100.0 / total << "%)" << boil::endl;
    boil::aout << "+------------------------------------------" << boil::endl;
  }
}

/***************************************************************************//**
*  Reports the time spent in the program with avereging over processors. 
*******************************************************************************/
void Timer :: report() {

  /* total time */
  real total = elapsed_s();
  boil::cart.average_real(& total);

  boil::oout << "+==========================================" << boil::endl;
  boil::oout << "| Total execution time: " << total << " [s]" << boil::endl;
  boil::oout << "+------------------------------------------" << boil::endl;

  if(total < 0.01) return;

  std::string elsew = "elsewhere";

  int longest(0);
  for(int i=0; i<sub.size(); i++) 
    longest=boil::maxi(longest, sub[i].name.length());
  for(int i=0; i<sub.size(); i++) 
    sub[i].name.resize(longest, ' ');
  elsew.resize(longest+3, ' '); /* 3 is length of "in " */

  real acum = 0.0; 
  for(int i=0; i<sub.size(); i++) {
    real here = sub[i].elapsed_s();
    boil::cart.average_real(& here);

    boil::oout << "| Time spent in " << sub[i].name << ": " 
               << here << " [s]    (" 
               << here * 100.0 / total << "%)" << boil::endl;
    acum += here;
  }
  if(sub.size()) {
    boil::oout << "| Time spent " << elsew << ": " 
               << total - acum << " [s]    ("
               << (total - acum) * 100.0 / total << "%)" << boil::endl;
    boil::oout << "+------------------------------------------" << boil::endl;
  }
}
