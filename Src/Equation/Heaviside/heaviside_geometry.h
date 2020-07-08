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
       real refval;
    };
    struct VAL3D {
      CELL3D cell;
      real value;
    };
    struct VERT {
      XYZ v;
      XYZ ref;
      real refval;
    };

    /* 2D */
    struct XY {
      real x,y;

      XY() {}
      XY(const real & a, const real & b) {
        x = a; y = b;
      }

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

      /* multiplication by number */
      XY operator *(const real & a) const {
        XY pv;
        pv.x = x*a;
        pv.y = y*a;

        return(pv);
      }

      /* division by number */
      XY operator /(const real & a) const {
        XY pv;
        pv.x = x/a;
        pv.y = y/a;

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
    
    /* the line object represent a plain line. However, for further
       operation (e.g. normal vector calculation), it is better to
       assume that it is actually an oriented line, i.e. vector from
       p[0] to p[1]. The methods in this class follow this assumption;
       moreover, convention is adopted that the normal vector to the 
       line points to the right of the line. */
    struct LINE {
      XY p[2];
      XY tangent; /* tangent vector */
      XY normal;  /* normal vector */
      real alpha; /* line constant */

      LINE(const XY & p1, const XY & p2) {
        p[0] = p1;
        p[1] = p2;

        assert(length()>0.);
        tangent.x = (p[1].x-p[0].x)/length();
        tangent.y = (p[1].y-p[0].y)/length();

        normal.x =  tangent.y;
        normal.y = -tangent.x;

        alpha = normal*p[0];
      }

      inline real length() const {
        return sqrt( (p[1].x-p[0].x)*(p[1].x-p[0].x)
                    +(p[1].y-p[0].y)*(p[1].y-p[0].y) );
      }

      inline XY CoM() const {
        return XY(0.5*(p[0].x+p[1].x),0.5*(p[0].y+p[1].y));
      }
    };
