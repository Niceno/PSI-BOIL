#ifndef FORMULA_H
#define FORMULA_H

#include "../Parallel/mpi_macros.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <cstring>

#include "scanner.h"
#include "../Parallel/Out/out.h"
#include "../Global/global_precision.h"

/***************************************************************************//**
*  \brief Evaluation of analytically prescribed expressions. 
*
*  Used for analytical prescription of boundary conditions, but its usage
*  might be envisaged for more general purposes: analytical prescription
*  of initial fields or maybe even some physical models. In later case,
*  the speed of interpretation might be a hinderance.
*
*  It's usage is probably best explained with two examples.
*
*  Example 1: Evaluation of a simple mathematical operation.
*  \code
*    f.evaluate("a=3");          // f.set("a",3); would do the same
*    f.evaluate("b=4");          // f.set("b",4); would do the same
*    boil::oout << f.evaluate("a+b") << boil::endl;     // prints 7
*  \endcode
*
*  Example 2: more elaborate example, initialize Scalar u with 
*             analytically prescribed function: u = x^2 + y^2 + z^2
*  \code
*    Domain d(gx, gx, gx);
*    Scalar u(d);
*  
*    // equation we would like to evaluate   
*    std::string equation("x*x + y*y + z*z");
*  
*    // browse through all the cells of variable "u"   
*    for_vijk(u, i,j,k) {
*  
*      // set equation variables   
*      f.set("x", u.xc(i));
*      f.set("y", u.yc(j));
*      f.set("z", u.zc(k));
*  
*      // evaluate equation   
*      u[i][j][k] = f.evaluate(equation);
*    }
*  \endcode
*
*******************************************************************************/

///////////////
//           //
//  Formula  //
//           //
///////////////
class Formula
 {
  public:
    // functions
    Formula()  {}
    ~Formula() {}
    real evaluate(std::string       & request);
    real evaluate(std::stringstream & request);
    real evaluate(const char * request);
    void set(const std::string &vname, real vvalue);
    void list();

  private:
    // data
    typedef struct var_struct
     {std::string label;
      real value;} var_type;

    char             token[256];
    int              ttype;
    std::vector<var_type> vars;
    const std::string emptyvar;
    char *sp;

    // functions
    real read_var(const std::string &vname);
    real expression();
    real simpleexpr();
    real term();
    real signedfactor();
    real function();
    real factor();
 };

#endif
