#ifndef ANGLE_NODE_H
#define ANGLE_NODE_H

/////////////////
//             //
//  AngleNode  // 
//             //
/////////////////
class AngleNode {

  public:
    AngleNode(const real an, const real xp, const real yp, const real zp);

    bool operator < (const AngleNode & other) const {
      return angle < other.angle;
    }

    friend std::ostream & 
      operator << (std::ostream & os, const AngleNode & an);

    real angle, x, y, z;
};

#endif
