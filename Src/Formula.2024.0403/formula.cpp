#include "formula.h"

/***************************************************************************//**
*  Evaluates string sent as a parameter.
*******************************************************************************/
real Formula :: evaluate(std::string & request)
/*-------------------+
|  calls expression  |
+-------------------*/
 {
  sp = &request[0];

  strcpy(token,"");
  if(0!=strpbrk(token,"+.")) ++sp; // what is this line doing ???
  scanner(&sp, token, ttype);
  return(expression());                         // <= expression()
 }

/***************************************************************************//**
*  Evaluates stringstring sent as a parameter. Just the interface to 
*  evaluate(string).
*******************************************************************************/
real Formula :: evaluate(std::stringstream & request)
 {
  std::string s = request.str();
  return evaluate(s);
 }

/***************************************************************************//**
*  Evaluates character array sent as a parameter.
*******************************************************************************/
real Formula :: evaluate(const char * request)
/*-------------------+
|  calls expression  |
+-------------------*/
 {
  std::string s(request);
  return evaluate(s);
 }

/******************************************************************************/
real Formula :: expression()
/*-------------------------+
|  calls simpleexpression  |
+-------------------------*/
 {
  real result=0;

  result=simpleexpr();
  while(ttype==PLUSOP || ttype==MINUSOP)
   {
    switch(ttype)
     {
      case PLUSOP : scanner(&sp, token, ttype);
   		    result+=simpleexpr();       // <= simpleexpression()
                    break;
      case MINUSOP: scanner(&sp, token, ttype);
   		    result-=simpleexpr();       // <= simpleexpression()
                    break;
     }
   }
  return(result);
 }

/******************************************************************************/
real Formula :: simpleexpr()
/*-------------+
|  calls term  |
+-------------*/
 {
  real result=0;

  result=term();
  while(ttype==MUXOP || ttype==DIVOP)
   {
    switch(ttype)
     {
      case MUXOP : scanner(&sp, token, ttype);
 		   result*=term();             // <= term()
                   break;
      case DIVOP : scanner(&sp, token, ttype);
 		   result/=term();             // <= term()
                   break;
     }
   }
  return(result);
 }

/******************************************************************************/
real Formula :: term(void)
/*---------------------+
|  calls signedfactor  |
+---------------------*/
 {
  real result=0;

  result=signedfactor();
  while(ttype==POWOP)
   {
    scanner(&sp, token, ttype);
    //result=exp(log(result)*signedfactor()); // <= signedfactor()
    /*---------------------------------------------------------------+ 
    |  If base is negative and exponent is not an integral value,    |
    |  or if base is zero and exponent is negative, a domain error   |
    |  occurs, setting the global variable errno to the value EDOM   |
    +---------------------------------------------------------------*/
    //result=pow( result, (int)signedfactor() ); // <= signedfactor()  
    result=pow( result, signedfactor() ); // <= signedfactor()  
   }
  return(result);
 }

/******************************************************************************/
real Formula :: signedfactor()
/*-----------------+
|  calls function  |
+-----------------*/
 {
  if(*token=='-')
   {
    scanner(&sp, token, ttype);
    return (0.0-function()); // <= function()
   }
  else if(*token=='+')
   {
    scanner(&sp, token, ttype);
    return function();       // <= function()
   }
  else
    return function();       // <= function()
 }

/******************************************************************************/
real Formula :: function() {
/*---------------+
|  calls factor  |
+---------------*/

  switch(ttype)
   {
    case SINF:  scanner(&sp, token, ttype);
                return sin(factor());      // <= factor()
    case COSF:  scanner(&sp, token, ttype);
                return cos(factor());      // <= factor()
    case TANF:  scanner(&sp, token, ttype);
                return tan(factor());      // <= factor()
    case ASINF: scanner(&sp, token, ttype);
                return asin(factor());     // <= factor()
    case ACOSF: scanner(&sp, token, ttype);
                return acos(factor());     // <= factor()
    case ATANF: scanner(&sp, token, ttype);
                return atan(factor());     // <= factor()
    case SQRT:  scanner(&sp, token, ttype);
                return sqrt(factor());     // <= factor()
   }

  return factor();
}


/******************************************************************************/
real Formula :: factor()
/*-------------------------------------------+
|  recursivelly calls: expression if needed  |
+-------------------------------------------*/
 {
  static real result=0;
  char holdvar[9];

  switch(ttype)
   {
    case REALLIT : result = atof(token); break;
    case EXPLIT  : result = atof(token); break;
    case INTLIT  : result = (real)atol(token); break;
    case LPAREN  : scanner(&sp, token, ttype);
                   result = expression(); break;
    case SYMBOL  : // pronasao je varijablu, staru ili novu
      strcpy(holdvar,token);
      scanner(&sp, token, ttype);

      if(ttype==ASSIGNOP) // ako se pridruzuje neka vrijednost varijabli
       {
        scanner(&sp, token, ttype);
        result=expression();
        set(holdvar,result);
        break;
       }
      else // ako se nista ne pridruzuje, nego se cita iz varijable
       {
        result=read_var(holdvar);
        // set(holdvar,result);
        // free(holdvar);
        return (result);
       }
   }
  scanner(&sp, token, ttype);
  return(result);
 }

/***************************************************************************//**
*  Prints out the names of all variables and their values.                   ,
*******************************************************************************/
void Formula :: list()
 {
  for(unsigned int i=0;i<vars.size();i++)
    boil::aout << vars[i].label << " = " << vars[i].value << boil::endl;
 }

/***************************************************************************//**
*  \param vname  - variable name
*  \param vvalue - value
*
*  Associates the value of variable "vname" with the value deined in "vvalue",
*  if the variable is defined. If it is not defined, then is created before
*  association with the value.
*******************************************************************************/
void Formula :: set(const std::string &vname, real vvalue)
 {
  for(unsigned int i=0;i<vars.size();i++) // da li postoji varijabla "vname"
    if( vars[i].label == vname ) // ako je nasao varijablu
     {vars[i].value=vvalue;
      return;}

  var_struct tmp;
  tmp.value = vvalue;
  tmp.label = vname;
  vars.push_back(tmp);

  return;
 }

/***************************************************************************//**
*  Returns the value of variable "vname", if it is defined. If it is not 
*  defined, prints an error message.
*******************************************************************************/
real Formula :: read_var(const std::string &vname)
 {
  real vvalue;

  // look if the variable is allready defined
  for(unsigned int i=0;i<vars.size();i++)
    if( vars[i].label == vname ) // it is here
     {vvalue=vars[i].value;
      return(vvalue);}

  // variable not defined (yet)
  boil::aout << "error in formula: variable \"" 
             << vname << "\" not found!" << boil::endl;

  return 0;
 }
