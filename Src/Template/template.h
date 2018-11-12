/***************************************************************************//**
*  \brief A template header file for PSI-Boil                 
*
*  This file was created as a part of the PSI-Boil Tutorial, Volume 2: 
*  Development Manual. It's purpose is to outline components a PSI-Boil
*  header file should have. It can be used as a starting point when defining
*  new classes.
*
*  \warning
*  Although it is compiled with the package, it has no functional purposes.
*******************************************************************************/

/* directives for the compiler */
#ifndef TEMPLATE_H   
#define TEMPLATE_H

#include <iostream>                     /* first stadard C++ headers  ...     */
#include "../Global/global_precision.h" /* ... followed by PSI-Boil's headers */

////////////////
//            //
//  Template  //
//            //
////////////////
class Template {
  public:
    //! Default construtor.
    Template() {at=0; bt=0;}  /* if short, it can stay in the header */
    
    //! A more elaborate constructor.
    Template(const real & x,
             const real & y); /* if long, defined in a .cpp source file */
    
    ~Template();

    //! Getter and setter functions should be in the header
    real a() const {return at;}    /* getter should be declared as const */ 
    real b() const {return bt;}    /* getter should be declared as const */ 
    void a(const real & x) {at=x;} /* const promoted for parameters */ 
    void b(const real & y) {bt=y;} /* const promoted for parameters */ 
                    
    //! Returns the sum of squares. 
    real sum_squares() const;

  private:
    real at, bt;        
};

#endif /* directive for the compiler */
