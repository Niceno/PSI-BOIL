#ifndef PLANE_H
#define PLANE_H

/////////////
//         //
//  Plane  //
//         //
/////////////
class Plane {

  public:
    Plane(const real xc, const real yc, const real zc, 
          const real nx, const real ny, const real nz) {
      assert( approx( (nx*nx + ny*ny + nz*nz), 1.0 ) );
      x0 = xc;
      y0 = yc;
      z0 = zc;
      A  = nx;
      B  = ny;
      C  = nz;
      D  = -A*x0 - B*y0 - C*z0;
    }

    real n(const int i) const {
      assert( i>=0 && i<= 2 );
      if(i==0) return A;
      if(i==1) return B;
      if(i==2) return C;
      return -1;
    }

    real x(const real yl, const real zl) const {
      if( !approx(A, 0.0) ) return (-B*yl - C*zl - D) / A;
      else                  return -FLT_MAX;
    }
    real y(const real xl, const real zl) const {
      if( !approx(B, 0.0) ) return (-A*xl - C*zl - D) / B;
      else                  return -FLT_MAX;
    }
    real z(const real xl, const real yl) const {
      if( !approx(C, 0.0) ) return (-A*xl - B*yl - D) / C;
      else                  return -FLT_MAX;
    }
    real distance(const real xx, const real yy, const real zz) const {
      return A*xx+B*yy+C*zz+D;
    }
  private:
    real x0;
    real y0;
    real z0;

    real A, B, C, D;
};

#endif
