/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (in Grid1D).
*
*  Used only inside Grid1D to denote the stride for multigrid coarsening.
*******************************************************************************/
#ifndef STEP_H
#define STEP_H

////////////
//        //
//  Step  //
//        //
////////////
class Step {
  public:
    explicit Step(const int i) : val(i) {}
    int size() const     {return val;}

  private:
    const int val;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: step.h,v 1.1 2009/07/17 16:35:20 niceno Exp $'/
+-----------------------------------------------------------------------------*/
