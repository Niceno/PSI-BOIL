#ifndef POLYGON_H
#define POLYGON_H

///////////////
//           //
//  Polygon  //
//           //
///////////////
class Polygon {

  public:
    Polygon(const int numb, 
            const real xp[], 
            const real yp[], 
            const real zp[], 
            const real np[] = NULL);

    //! Number of polygon's nodes (three for triangle, for example).
    int nnodes() const {return nn;}
 
    real xn(const int i) const {return x  [i];}
    real yn(const int i) const {return y  [i];}
    real zn(const int i) const {return z  [i];}
    real n (const int i) const {return nor[i];}

    //! Center of gravity.
    real xg() const {return g[0];}
    real yg() const {return g[1];}
    real zg() const {return g[2];}

    real area_x() const;
    real area_y() const;
    real area_z() const;

    real a() const {return A;}
    real b() const {return B;}
    real c() const {return C;}
    real d() const {return D;}

    real minx() const {return min[0];}
    real maxx() const {return max[0];}
    real miny() const {return min[1];}
    real maxy() const {return max[1];}
    real minz() const {return min[2];}
    real maxz() const {return max[2];}

    bool foot(const real a, const real b, const real c, real & d);

    void split(Polygon ** left, Polygon ** right) const;
    void arrange();

  private:
    int  nn;   /* number of nodes */
    real x[6];
    real y[6];
    real z[6];
    real nor[3];
    real g[3];   /* center of gravity (still not used) */
    real max[3]; /* max x, y and z cooridnates */
    real min[3]; /* min x, y and z cooridnates */

    real A, B, C, D;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: polygon.h,v 1.3 2015/08/17 10:37:05 niceno Exp $'/
+-----------------------------------------------------------------------------*/
