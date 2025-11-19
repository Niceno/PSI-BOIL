///////////////////////
//                   //
//  Cycle selection  //
//                   //
///////////////////////
/* this is a ravioli class for multigrid cycle selection */
class Cycle {
  public:
    Cycle() {val=-1;}

    static const Cycle none() {return Cycle(-1);}
    static const Cycle Z()    {return Cycle(-1);}
    static const Cycle V()    {return Cycle( 1);}
    static const Cycle F()    {return Cycle( 2);}
    static const Cycle W()    {return Cycle( 4);}
    static const Cycle flex() {return Cycle( 5);}

    //! Prints the cycle name.
    friend std::ostream & operator << (std::ostream & ost, const Cycle & com) {
      switch(com.val) {
        case(-1): ost << "none"; break;
        case( 1): ost << "V-cycle"; break;
        case( 2): ost << "F-cycle"; break;
        case( 4): ost << "W-cycle"; break;
        case( 5): ost << "flex-cycle"; break;
      }

      return ost;
    }

    bool operator == (const Cycle & o) const {return val == o.val;}
    bool operator != (const Cycle & o) const {return val != o.val;}

  private:
    int val;

    /* avoid implicit conversions of integer to Cycle */
    explicit Cycle(const int m) {val = m;}
};
