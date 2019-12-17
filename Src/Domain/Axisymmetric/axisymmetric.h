#ifndef AXISYMMETRIC_H
#define AXISYMMETRIC_H

#include "../domain.h"

/***************************************************************************//**
* A derived class from the Domain class  
* - x-direction is the radial direction
* - z-direction is the axial direction
* - y-direction is a dummy direction
*******************************************************************************/
///////////////////////////
//                       //
//  Axisymmetric domain  //
//                       //
///////////////////////////
class Axisymmetric : public Domain {
  public:
    Axisymmetric(const Grid1D & ogx, const Grid1D & ogz, const real dy,
                 Body * b = NULL, /* it will change, that is why it is pointer */
                 const std::string n="axisym_domain", 
                 const Decompose dec=Decompose::xz(),
                 const bool print_statistics = true) :
    Domain(ogx,Grid1D(dy),ogz,b,n,dec,false),
    ydummy(dy),
    angle(1.0) { 
      grid_y_original = &ydummy;
      check_radial_grid(ogx);

      if(b) {
        if(print_statistics)
          statistics(body);
      } else {
        if(print_statistics)
          statistics();
      }

    };

    ~Axisymmetric() {};

    /* get angle of simulated wedge in radians */
    real wedge_angle() { return angle; }

    /* computes and prints grid statistics
     * -> virtual so that it gets called from the proper constructor */
    virtual void statistics(Body * b = NULL);

    /* checks */
    virtual bool is_axisymmetric() const { return true; }
    virtual bool is_cartesian() const { return false; }


    /* cell surfaces */
    virtual real dSx(const int i, const int j, const int k) const;
    virtual real dSy(const int i, const int j, const int k) const;
    virtual real dSz(const int i, const int j, const int k) const;

    virtual real dSx_xstag(const int i, const int j, const int k) const;
    virtual real dSx_ystag(const int i, const int j, const int k) const;
    virtual real dSx_zstag(const int i, const int j, const int k) const;

    virtual real dSy_xstag(const int i, const int j, const int k) const;
    virtual real dSy_ystag(const int i, const int j, const int k) const;
    virtual real dSy_zstag(const int i, const int j, const int k) const;

    virtual real dSz_xstag(const int i, const int j, const int k) const;
    virtual real dSz_ystag(const int i, const int j, const int k) const;
    virtual real dSz_zstag(const int i, const int j, const int k) const;

    virtual real dSx(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSy(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSz(const Sign sig, const int i, const int j, const int k) const;

    virtual real dSx_xstag(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSx_ystag(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSx_zstag(const Sign sig, const int i, const int j, const int k) const;

    virtual real dSy_xstag(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSy_ystag(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSy_zstag(const Sign sig, const int i, const int j, const int k) const;

    virtual real dSz_xstag(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSz_ystag(const Sign sig, const int i, const int j, const int k) const;
    virtual real dSz_zstag(const Sign sig, const int i, const int j, const int k) const;

    /* cell volume */
    virtual real dV(const int i, const int j, const int k) const;

    virtual real dV_xstag(const int i, const int j, const int k) const;
    virtual real dV_ystag(const int i, const int j, const int k) const;
    virtual real dV_zstag(const int i, const int j, const int k) const;

  private:
    void check_radial_grid(const Grid1D & gx);

    Grid1D ydummy; 	  
    real angle; /* angle of wedge. Do not change without modifying geometrical functions!*/
};

#endif
