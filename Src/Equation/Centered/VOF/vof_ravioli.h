////////////////////////////
//                        //
//  Normal vector method  //
//                        //
////////////////////////////
/* this is a ravioli class for normal vector method selection */
class NormMethod {
  public:
    NormMethod() {val=-1;}

    static const NormMethod undefined() {return NormMethod(-1);}
    static const NormMethod Young()     {return NormMethod( 1);}
    static const NormMethod Mixed()     {return NormMethod( 2);}
    static const NormMethod CC()        {return NormMethod( 3);}
    static const NormMethod ElviraXZ()  {return NormMethod( 4);}
    static const NormMethod ElviraXY()  {return NormMethod( 5);}
    static const NormMethod ElviraYZ()  {return NormMethod( 6);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const NormMethod & com) {
      switch(com.val) {
        case(-1): ost << "undefined"; break;
        case( 1): ost << "Young"; break;
        case( 2): ost << "Mixed"; break;
        case( 3): ost << "CC"; break;
        case( 4): ost << "ElviraXZ"; break;
        case( 5): ost << "ElviraXY"; break;
        case( 6): ost << "ElviraYZ"; break;
      }

      return ost;
    }

    bool operator == (const NormMethod & o) const {return val == o.val;}
    bool operator != (const NormMethod & o) const {return val != o.val;}

  private:
    int val;

    /* avoid implicit conversions of integer to Norm */
    explicit NormMethod(const int m) {val = m;}
};

///////////////////
//               //
//  Curv method  //
//               //
///////////////////
/* this is a ravioli class for curvature method selection */
class CurvMethod {
  public:
    CurvMethod() {val=-1;}

    static const CurvMethod undefined() {return CurvMethod(-1);}
    static const CurvMethod DivNorm()   {return CurvMethod( 1);}
    static const CurvMethod HF()        {return CurvMethod( 2);}
    /* 2D method, X = wall tangent, Z = wall normal dir */
    static const CurvMethod HFmixedXZ() {return CurvMethod( 3);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const CurvMethod & com) {
      switch(com.val) {
        case(-1): ost << "undefined"; break;
        case( 1): ost << "DivNorm"; break;
        case( 2): ost << "HF"; break;
        case( 3): ost << "HFmixedXZ"; break;
      }

      return ost;
    }

    bool operator == (const CurvMethod & o) const {return val == o.val;}
    bool operator != (const CurvMethod & o) const {return val != o.val;}

  private:
    int val;

    /* avoid implicit conversions of integer to Curv */
    explicit CurvMethod(const int m) {val = m;}
};

////////////////////////
//                    //
//  Detachment model  //
//                    //
////////////////////////
/* this is a ravioli class for detachment modeling */
class DetachmentModel {
  public:
    DetachmentModel() :
      init(false),
      detach(false),
      etacrit(boil::unreal),
      etacurr(boil::unreal)
    {}

    DetachmentModel(const real cangle) : 
      init(true),
      detach(false),
      etacurr(boil::unreal)
    { 
      set_detachment_params(cangle);
    }

    inline bool initialized() const {
      return init;
    }

    inline bool detached() const {
      return detach;
    }

    inline void test_detachment(const real eta) {
      etacurr = eta;
      bool test = (eta<=etacrit);
      if(detach != test) {
        boil::oout<<"Detachment status changed! Was/is: "
                  <<detach<<" "<<test<<" | "
                  <<etacurr<<" < "<<etacrit<<boil::endl;
        detach = test;
      }
      return;
    }

    void set_detachment_params(const real cangle) {
      init = true;
      cangleval = cangle;
      real cval = std::max(boil::atto,std::min(cangle,90.));
      /* dimensionless detachment height */
      etacrit = 0.5*cos(cval*boil::pi/180.)/sin(cval*boil::pi/180.);

      boil::oout<<"VOF: detachment model invoked with critical height= "
                 <<etacrit<<" (angle "<<cval<<")"
                 <<boil::endl;

    }

  private:
    bool init, detach;
    real cangleval;
    real etacrit, etacurr;
};