#ifndef TWOLEVEL_H
#define TWOLEVEL_H

#include "../Axisymmetric/axisymmetric.h"

////////////////////////
//                    //
//  Two-level domain  //
//                    //
////////////////////////
class TwoLevelDomain {
  public:
    Domain & coarse() { return *cd; }
    Domain & fine()   { return *fd; }
    const Domain & coarse() const { return *cd; }
    const Domain & fine() const   { return *fd; }

  protected:
    TwoLevelDomain(Body * b = NULL)  
    {
      if(b) {
        cbody = new Body(b->name);
      } else {
        cbody = NULL;
      }
    }

    ~TwoLevelDomain() { delete cd; delete fd; delete cbody; }

    Body * cbody;
    Domain * cd, * fd;

};

////////////////////////
//                    //
//  Two-level domain  //
//                    //
////////////////////////
class TwoLevelCartesian : public TwoLevelDomain {
  public:
    TwoLevelCartesian(const Grid1D & ogx, const Grid1D & ogy, const Grid1D & ogz,
                   Body * b = NULL, 
                   const std::string n="domain",
                   const DimCut dmc = DimCut::xyz, 
                   const Decompose dec = Decompose::xyz(),
                   const bool print_statistics = true) :
      TwoLevelDomain(b)
    {

      fd = new Domain(ogx,ogy,ogz,b,n,dec,print_statistics);

      Step cx(2), cy(2), cz(2);

      switch(dmc) {
        case DimCut::x   : cy = Step(1); cz = Step(1); break;
        case DimCut::y   : cx = Step(1); cz = Step(1); break;
        case DimCut::z   : cx = Step(1); cy = Step(1); break;
        case DimCut::xy  : cz = Step(1); break;
        case DimCut::xz  : cy = Step(1); break;
        case DimCut::yz  : cx = Step(1); break;
        default : break;
      }
      cd = new Domain(*fd,cx,cy,cz,cbody,print_statistics);
    }
};

/////////////////////////////////////
//                                 //
//  Two-level axisymmetric domain  //
//                                 //
/////////////////////////////////////
class TwoLevelAxisymmetric : public TwoLevelDomain {
  public:
    TwoLevelAxisymmetric(const Grid1D & ogx, const Grid1D & ogz, const real dy,
                         Body * b = NULL,
                         const std::string n="axisym_domain",
                         const DimCut dmc = DimCut::xz,
                         const Decompose dec = Decompose::xz(),
                         const bool print_statistics = true) :
      TwoLevelDomain(b)
    {
      fd = new Axisymmetric(ogx,ogz,dy,b,n,dec,print_statistics);

      Step cx(2), cz(2);

      switch(dmc) {
        case DimCut::x : cz = Step(1); break;
        case DimCut::z : cx = Step(1); break;
        default : break;
      }
      cd = new Axisymmetric(*fd,cx,cz,cbody,print_statistics);
    }
};
#endif
