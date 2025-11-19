#ifndef SIMULATION_TIME_H
#define SIMULATION_TIME_H

#define CFL_HIST_SIZE 16

#include <cfloat>

#include "../Global/global_precision.h"
#include "../Global/global_malloc.h"
#include "../Parallel/Out/print.h"
#include "../Timer/timer.h"

/////////////
//         //
//  Times  //
//         //
/////////////
class Times {
  public:
    /* nothing is specified -> steady */
    explicit Times();

    /* number of time steps and time step are specified */
    explicit Times(const int n, const real d);

    /* total time and time step are specified */
    explicit Times(const real t, const real d);

    /* total time and number of time steps are specified */
    explicit Times(const real t, const int n);

    int  total_steps()  const {return t_dt;}
    int  current_step() const {return c_dt + f_dt;}
    int  first_step()   const {return f_dt;}
    int  first_step(int f)    {return f_dt = f;}
    real dt()           const {return d_t;}
    real dti()          const {return d_ti;}
    real total_time()   const {return t_time;}
    real current_time() const {return time;}
    void current_time(const real t) {time = t;}
    void start()    {c_dt = 1;}
    bool end()      {
                     if((time<=t_time) && ((c_dt+f_dt)<=t_dt))  
                      {print_current(); 
                       return true;}
                     else
                       return false;
                    } 
    void increase() {
                       c_dt++;
                       real tnew = 0.0;
                       if(boil::cart.iam()==0){
                         tnew=time+d_t;
                       }
                       boil::cart.sum_real(&tnew);
                       time = tnew;
                     }
    void set_dt(const real d) {d_t = d; d_ti=1.0/d_t;}
    void increase_nhist() {nhist++;}
    void decrease_nhist() {nhist--;}
    int  get_nhist() {return nhist;}
    real get_cfl_hist(int i) {return cfl_hist[i];}
    void set_cfl_hist(int i, real r) {cfl_hist[i]=r;}
    void control_dt(const real current_cfl, 
                    const real target_cfl, 
                    const real max_dt = FLT_MAX);
    void set_coef_inc(real r){boil::oout<<"simulation_time:coef_inc= "<<r<<"\n";
                              coef_inc=r;}
    void set_coef_dec(real r){boil::oout<<"simulation_time:coef_dec= "<<r<<"\n";
                              coef_dec=r;}
    void print_time(bool b) {print_cinfo=b;}

  private:
    void print_total() const;
    void print_current() const;
    int  f_dt;   /* fist time step (used in restart) */
    int  c_dt;   /* current time step in this run */
    real time;   /* current time for this run */
    int  t_dt;   /* total number of time steps in this run */
    real d_t;    /* time step */
    real d_ti;   /* inverse of the time step */
    real t_time; /* total time */
    real * cfl_hist; /* history of cfl */
    real coef_inc, coef_dec; /* coefficient for control_dt */
    int nhist, istep_inc;
    bool print_cinfo;
};

#endif
