#include "floodfill.h"

/***************************************************************************//**
*  scans the face of the decomposed domain storing values on either side
*  of the face
*******************************************************************************/
void Floodfill::face_scan_rgn() {

  int in_rgns, out_rgns;

  /* scan both faces in i direction */
  for_vjk(rgnid,j,k){
    in_rgns  = int(rgnid[rgnid.ei()][j][k]);
    out_rgns = int(rgnid[rgnid.ei()+1][j][k]);
    compare_rgn(in_rgns, out_rgns);
  }
  for_vjk(rgnid,j,k){
    in_rgns  = int(rgnid[rgnid.si()][j][k]);
    out_rgns = int(rgnid[rgnid.si()-1][j][k]);
    compare_rgn(in_rgns, out_rgns);
  }

  /* scan both faces in j direction */
  for_vik(rgnid,i,k){
    in_rgns  = int(rgnid[i][rgnid.ej()][k]);
    out_rgns = int(rgnid[i][rgnid.ej()+1][k]);
    compare_rgn(in_rgns, out_rgns);
  }
  for_vik(rgnid,i,k){
    in_rgns  = int(rgnid[i][rgnid.sj()][k]);
    out_rgns = int(rgnid[i][rgnid.sj()-1][k]);
    compare_rgn(in_rgns, out_rgns);
  }

  /* scan both faces in k direction */
  for_vij(rgnid,i,j){
    in_rgns  = int(rgnid[i][j][rgnid.ek()]);
    out_rgns = int(rgnid[i][j][rgnid.ek()+1]);
    compare_rgn(in_rgns, out_rgns);
  }
  for_vij(rgnid,i,j){
    in_rgns  = int(rgnid[i][j][rgnid.sk()]);
    out_rgns = int(rgnid[i][j][rgnid.sk()-1]);
    compare_rgn(in_rgns, out_rgns);
  }
}


/*-----------------------------------------------------------------------------+
 '$Id: identify_regions  lafferty
+-----------------------------------------------------------------------------*/
