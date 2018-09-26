#include "body.h"

/******************************************************************************/
real prod_inner(const real vec1[], const real vec2[]) {
  return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}

/******************************************************************************/
real vnorm(const real vec[]) {
  return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

/******************************************************************************/
void prod_outer(const real vec1[], const real vec2[], real out[]) {
  out[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  out[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  out[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  return;
}

/******************************************************************************/
bool Polygon::foot(const real xin, const real yin, const real zin,
                   real & alen) {
/***************************************************************************//**
* /brief calculate foot of perpendicular.
*        return true if foot is inside the polygon
*        return false if foot is outside the polygon
*        output : alen (length of perpendicular line)
*    argorithm for determination of inside/outside
*      set the foot at orign.
*      project edge of polygon to the circle of radius = 1
*      sum ((angle of view) * sign )
*      if the foot is outside, summantion is zero.
*******************************************************************************/
      
  /* coordinates of foot */
  real t = -(A*xin+B*yin+C*zin+D)/(A*A+B*B+C*C);
  real xft = xin + t*A;
  real yft = yin + t*B;
  real zft = zin + t*C;

  /* distance of perpendicular line */
  alen = sqrt ( pow(xin-xft,2.0)
              + pow(yin-yft,2.0)
              + pow(zin-zft,2.0));

  /* check foot is inside/outside of Polygon */
  real sum_angle=0.0;
  real vec1[3],vec2[3];
  real norm_m[3];
  for(int iv=0; iv<nn-1; iv++){
    vec1[0]=xn(iv)-xft;
    vec1[1]=yn(iv)-yft;
    vec1[2]=zn(iv)-zft;
    vec2[0]=xn(iv+1)-xft;
    vec2[1]=yn(iv+1)-yft;
    vec2[2]=zn(iv+1)-zft;
    if(vnorm(vec1)==0||vnorm(vec2)==0) return true;
    real ac=std::max(-1.0,prod_inner(vec1,vec2)/(vnorm(vec1)*vnorm(vec2)));
    ac=std::min(ac,1.0);
    real theta = acos(ac);
    //avoid truncation error: intel compiler faces floating point exception.
    //real theta = acos(prod_inner(vec1,vec2)/(vnorm(vec1)*vnorm(vec2)));
    prod_outer(vec1,vec2,norm_m);
    real asign = copysign(1.0,norm_m[0]*A+norm_m[1]*B+norm_m[2]*C);
    sum_angle += asign*theta;
  }
    vec1[0]=xn(nn-1)-xft;
    vec1[1]=yn(nn-1)-yft;
    vec1[2]=zn(nn-1)-zft;
    vec2[0]=xn(0)-xft;
    vec2[1]=yn(0)-yft;
    vec2[2]=zn(0)-zft;
    if(vnorm(vec1)==0||vnorm(vec2)==0) return true;
    real ac=std::max(-1.0,prod_inner(vec1,vec2)/(vnorm(vec1)*vnorm(vec2)));
    ac=std::min(ac,1.0);
    real theta = acos(ac);
    //avoid truncation error: intel compiler faces floating point exception.
    //real theta = acos(prod_inner(vec1,vec2)/(vnorm(vec1)*vnorm(vec2)));
    prod_outer(vec1,vec2,norm_m);
    real asign = copysign(1.0,norm_m[0]*A+norm_m[1]*B+norm_m[2]*C);
    sum_angle += asign*theta;

  if(fabs(sum_angle)>2.0*3.1415){  //2*pi
    return true;
  } else {
    return false;
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: polygon_foot.cpp,v 1.2 2013/10/02 15:56:54 sato Exp $'/
+-----------------------------------------------------------------------------*/
