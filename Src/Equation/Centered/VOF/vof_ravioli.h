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
    static const CurvMethod none()      {return CurvMethod( 0);}
    static const CurvMethod DivNorm()   {return CurvMethod( 1);}
    static const CurvMethod HF()        {return CurvMethod( 2);}
    /* 2D methods, X = wall tangent, Z = wall normal dir */
    static const CurvMethod HFparallelXZ() {return CurvMethod( 3);}
    static const CurvMethod HFnormalXZ()   {return CurvMethod( 4);}
    static const CurvMethod HFmixedXZ()    {return CurvMethod( 5);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const CurvMethod & com) {
      switch(com.val) {
        case(-1): ost << "undefined"; break;
        case( 0): ost << "none"; break;
        case( 1): ost << "DivNorm"; break;
        case( 2): ost << "HF"; break;
        case( 3): ost << "HFparallelXZ"; break;
        case( 4): ost << "HFnormalXZ"; break;
        case( 5): ost << "HFmixedXZ"; break;
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

///////////////////////
//                   //
//  Topology method  //
//                   //
///////////////////////
/* this is a ravioli class for topology method selection */
class TopoMethod {
  public:
    TopoMethod() {val=-1;}

    static const TopoMethod undefined() {return TopoMethod(-1);}
    static const TopoMethod Hybrid()    {return TopoMethod( 1);}
    static const TopoMethod Heaviside() {return TopoMethod( 2);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const TopoMethod & com) {
      switch(com.val) {
        case(-1): ost << "undefined"; break;
        case( 1): ost << "Hybrid"; break;
        case( 2): ost << "Heaviside"; break;
      }

      return ost;
    }

    bool operator == (const TopoMethod & o) const {return val == o.val;}
    bool operator != (const TopoMethod & o) const {return val != o.val;}

  private:
    int val;

    /* avoid implicit conversions of integer to Topo */
    explicit TopoMethod(const int m) {val = m;}
};

//////////////////////
//                  //
//  Subgrid method  //
//                  //
//////////////////////
/* this is a ravioli class for subgrid method selection */
class SubgridMethod {
  public:
    SubgridMethod() {val=-1;}

    static const SubgridMethod undefined()  {return SubgridMethod(-1);}
    static const SubgridMethod None()       {return SubgridMethod( 1);}
    static const SubgridMethod PLIC()       {return SubgridMethod( 2);}
    static const SubgridMethod SLICliquid() {return SubgridMethod( 3);}
    static const SubgridMethod SLICgas()    {return SubgridMethod( 4);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const SubgridMethod & com) {
      switch(com.val) {
        case(-1): ost << "undefined"; break;
        case( 1): ost << "None"; break;
        case( 2): ost << "PLIC"; break;
        case( 3): ost << "SLIC-liquid"; break;
        case( 4): ost << "SLIC-gas"; break;
      }

      return ost;
    }

    bool operator == (const SubgridMethod & o) const {return val == o.val;}
    bool operator != (const SubgridMethod & o) const {return val != o.val;}

  private:
    int val;

    /* avoid implicit conversions of integer to Subgrid */
    explicit SubgridMethod(const int m) {val = m;}
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

/////////////
//         //
//  HFset  //
//         //
/////////////
/* this is a ravioli struct for holding HF parameters,
 *
 *   mof, nof, majorext, minorext --- (symmetric) grid size
 *   iterloop --- extrapolation
 *   theta_crit --- for 3D, Eq. (6) in Lopez's paper
 *   blending_angle, n0square, n0, overshoot --- for 2D blending, 
     unused at the moment (2020.02.28) */
struct HFset {
  HFset() : mof(3), nof(1), /* symmetric stencil is constructed */
            majorext(mof*2+1), minorext(nof*2+1),
            /* warning: if stencil size is changed in the minorext direction,
                        the kernel must be properly adjusted! */
            iterloop(5), /* 2019.07.09 */
            theta_crit(0.8),
            blending_angle(40./180.*boil::pi),
            n0square(1.-1./(1.+tan(blending_angle)*tan(blending_angle))),
            n0(sqrt(n0square)),
            overshoot(0.5)
  {};

  const int mof, nof, majorext, minorext;
  const int iterloop;
  const real theta_crit;
  const real blending_angle, n0square, n0, overshoot;
};

