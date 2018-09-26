#include "region.h"

/******************************************************************************/
/* constructor */
Region::Region(const int rid,
               const int cvol,
               const real x, const real y, const real z,
               const real u, const real v, const real w) {

  m_id=rid;
  m_pos[0]=x; m_pos[1]=y; m_pos[2]=z;
  m_opos[0]=x; m_opos[1]=y; m_opos[2]=z;
  m_vel[0]=u; m_vel[1]=v; m_vel[2]=w;
  m_comvel[0]=u; m_comvel[1]=v; m_comvel[2]=w;
  m_cellvol=cvol;
  m_hiding=false;
  m_tsteps_hiding=0;

}

/*-----------------------------------------------------------------------------+
 '$Id: region.cpp,v 1.1 2018/02/16 19:07:11 sato Exp $'/ 
+-----------------------------------------------------------------------------*/
