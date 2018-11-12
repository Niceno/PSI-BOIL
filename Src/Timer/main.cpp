#include <cmath>
#include "timer.h"
#include "../Global/global_endl.h"

void fa();
void fb();
void fc();
void fe();

/******************************************************************************/
main(int argc, char * argv[]) {CART.init(argc, argv);

  timer.start();

  for(int i=0; i<3; i++) {
    fc();
    fe();
  }

  timer.stop();
  timer.report_separate();
  timer.report();
}

/******************************************************************************/
void fa() {
  real a,b=12.0;

  timer.start("a");

  for(int i=0; i<1000; i++) 
    for(int j=0; j<1000; j++)
      for(int k=0; k<1000; k++)
        a = sqrt(b);

  fb();

  timer.stop("a");
}

/******************************************************************************/
void fb() {
  real b=12.0;

  timer.start("b");

  for(int i=0; i<1000; i++) 
    for(int j=0; j<1000; j++)
      for(int k=0; k<1000; k++)
        real a = sqrt(b);

  timer.stop("b");
}

/******************************************************************************/
void fc() {
  timer.start("c");

  fa();

  timer.stop("c");
}

/******************************************************************************/
void fe() {
  real b=12.0;
  for(int i=0; i< 500; i++) 
    for(int j=0; j< 500; j++)
      for(int k=0; k< 500; k++)
        real a = sqrt(b);
}
