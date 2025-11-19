#include "enthalpytif.h"

/******************************************************************************/
real EnthalpyTIF::temperature_node(real len_s, real lam_s, real tmp_s
                                  , real len_f, real lam_f, real tmp_f) {
/***************************************************************************//**
*  \brief calculate temperature at node point
*             len_s         len_f
*             lam_s         lam_f
*         *-------------*------------*
*       tmp_s       tmp_node        tmp_f
*******************************************************************************/
  return (len_f*lam_s*tmp_s + len_s*lam_f*tmp_f)/(len_f*lam_s + len_s*lam_f);
}


