#ifndef HEAVISIDE_H
#define HEAVISIDE_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Domain/domain.h"

/////////////////
//             //
//  Heaviside  //
//             //
/////////////////
class Heaviside { /* this class is an abstract class! */
  public:
    Heaviside(const Scalar * CLR, Scalar * PHI = NULL, Scalar * ADENS = NULL,
              const real CLRSURF = 0.5) : 
      clr(CLR), dom((*CLR).domain()), phi(PHI), adens(ADENS), clrsurf(CLRSURF) {}
    ~Heaviside() {};

    const Domain * domain() const {return dom;}
    void calculate(const bool evalflag = true);
    void calculate_heaviside(const bool evalflag = true);
    void calculate_adens(const bool evalflag = true);

    /* pure virtual functions */
    virtual void evaluate_nodes() = 0;
    virtual real area(const int i, const int j, const int k) = 0;
    virtual real vf(const int i, const int j, const int k) = 0;

    real operator() (const int i, const int j, const int k) const {
      return (*phi)[i][j][k];
    }

    Scalar nodalvals;

  protected:
    const Scalar * clr;
    const Domain * dom; 
    Scalar * phi;
    Scalar * adens;
    real clrsurf;

    /* 3D */
    struct XYZ {
      real x,y,z;

      /* negation */
      XYZ operator -() const {
        XYZ nv;
        nv.x = -x;
        nv.y = -y;
        nv.z = -z;

        return(nv);
      }

      /* addition and subtraction */
      XYZ operator +(const XYZ & p2) const {
        XYZ pv;
        pv.x = x + p2.x;
        pv.y = y + p2.y;
        pv.z = z + p2.z;

        return(pv);
      }

      XYZ operator -(const XYZ & p2) const {
        XYZ pv;
        pv.x = x - p2.x;
        pv.y = y - p2.y;
        pv.z = z - p2.z;

        return(pv);
      }

      /* dot product */
      inline real operator *(const XYZ & p2) const {
        return x*p2.x + y*p2.y + z*p2.z; 
      }

      XYZ CrossProduct(const XYZ & p2) const {
        XYZ cp;
        cp.x = y * p2.z - z * p2.y;
        cp.y = z * p2.x - x * p2.z;
        cp.z = x * p2.y - y * p2.x;

        return(cp);
      }
    }; 

    struct TRIANGLE {
       XYZ p[3];
       XYZ v[3];

       real area() const {
         real x1 = p[1].x-p[0].x;
         real y1 = p[1].y-p[0].y;
         real z1 = p[1].z-p[0].z;
         real x2 = p[2].x-p[0].x;
         real y2 = p[2].y-p[0].y;
         real z2 = p[2].z-p[0].z;
         real val = (y1*z2-z1*y2)*(y1*z2-z1*y2)
                    +(z1*x2-x1*z2)*(z1*x2-x1*z2)
                    +(x1*y2-y1*x2)*(x1*y2-y1*x2);
         return 0.5 * sqrt(val);
       }

    };
    struct CELL3D {
       XYZ p[8];
       real val[8];
    };
    struct VAL3D {
      CELL3D cell;
      real value;
      real refval;
    };
    struct VERT {
      XYZ v;
      XYZ ref;
      real refval;
    };

    /* 2D */
    struct XY {
      real x,y;

      /* negation */
      XY operator -() const {
        XY nv;
        nv.x = -x;
        nv.y = -y;

        return(nv);
      }

      /* addition and subtraction */
      XY operator +(const XY & p2) const {
        XY pv;
        pv.x = x + p2.x;
        pv.y = y + p2.y;

        return(pv);
      }

      XY operator -(const XY & p2) const {
        XY pv;
        pv.x = x - p2.x;
        pv.y = y - p2.y;

        return(pv);
      }

      /* dot product */
      inline real operator *(const XY & p2) const {
        return x*p2.x + y*p2.y; 
      }
    };

    struct CELL2D {
      XY p[4];
      real val[4];
      real refval;
    };
    struct VAL2D {
      CELL2D cell;
      real value;
    };
    struct LINE {
      XY p[2];

      LINE(const XY & p1, const XY & p2) {
        p[0] = p1;
        p[1] = p2;
      }
    };

    /* interpolation */
    XY  VertexInterp(const real & isolevel,
                     const XY & p1, const XY & p2, 
                     const real & valp1, const real & valp2);
    XYZ VertexInterp(const real & isolevel,
                     const XYZ & p1, const XYZ & p2,
                     const real & valp1, const real & valp2);

    /* cross product */
    XYZ CrossProduct(const XYZ & p1, const XYZ & p2);

    /* shoelace */
    real Shoelace(const XY & v1, const XY & v2, const XY & v3);
    real Shoelace(const XY & v1, const XY & v2, const XY & v3,
                  const XY & v4);
    real Shoelace(const XY & v1, const XY & v2, const XY & v3,
                  const XY & v4, const XY & v5);

    /* surface divergence */ 
    real triangle_surface_divergence(const TRIANGLE & t);

    /* main body of marching squares*/
    real standing_square(const CELL2D & grid, const real & isolevel,
                         const real & totarea, std::vector<LINE> & lines);
};

#endif
