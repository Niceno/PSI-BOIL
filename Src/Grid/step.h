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

    //! Prints the step value.
    friend std::ostream & operator << (std::ostream & ost, const Step & stp) {
      ost << stp.size();
      return ost;
    }

    const Step & operator = (const Step & s)
     {val = s.size(); return *this;}

  private:
    int val;
};

#endif
