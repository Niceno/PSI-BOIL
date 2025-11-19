#include "floodfill.h"

/***************************************************************************//**
*  Assigns closed regions unique integer identifier number in Scalar rgn   
*******************************************************************************/

void Floodfill::scanline() {

  int rgn_inc;
  int i=new_seed_stack.back().i();
  int j=new_seed_stack.back().j();
  int k=new_seed_stack.back().k();
  new_seed_stack.pop_back(); 
  int j1;

  bool spanw, spane, spanb, spant;
  if (!same_fill_stack.empty()) {
    boil::aout<<"EXIT: same_fill_stack not empty in FF_scanline()" 
      << boil::endl;
    exit(1);
  }
  if (rgnid[i][j][k] != 0) return;  
  same_fill_stack.push_back(Index_ijk(i, j, k));
  if (c[i][j][k] <= c_criterion) {
    m_rgn_n--;
    rgn_inc = m_rgn_n;
  }
  else if (c[i][j][k] > c_criterion) {
    m_rgn_p++;
    rgn_inc = m_rgn_p;
  }
  while (!same_fill_stack.empty()) {
    Index_ijk same_fill_stackpt = same_fill_stack.back();
    same_fill_stack.pop_back();  //remove last element from vector stack
    j1 = same_fill_stackpt.j();
    i = same_fill_stackpt.i();
    k = same_fill_stackpt.k();
  
    while ((j1 >= c.sj()) && (rgnid[i][j1][k]==0)) {
      if (((rgn_inc < 0) && (c[i][j1][k] <= c_criterion)) ||
          ((rgn_inc > 0) && (c[i][j1][k] > c_criterion))) {
        j1--;
      }
      else {
        new_seed_stack.push_back(Index_ijk(i,j1,k));
        break;
      }
    }
    j1++;
  
    spanw = spane = spanb = spant = 0;
    while ((j1 <= c.ej()) && (rgnid[i][j1][k]==0)) {
      if (((rgn_inc < 0) && (c[i][j1][k] <= c_criterion)) ||  
          ((rgn_inc > 0) && (c[i][j1][k] > c_criterion))) {
        rgnid[i][j1][k] = rgn_inc;  //rgnid[i][j][k] assigned as rgn_inc
        if ((!spanw && i>c.si()) && (rgnid[i-1][j1][k]==0)) {
          if (((rgn_inc < 0) && (c[i-1][j1][k] <= c_criterion)) ||  
              ((rgn_inc > 0) && (c[i-1][j1][k] > c_criterion))) {
            
            same_fill_stack.push_back(Index_ijk(i-1,j1,k));
            spanw = 1;
          }
          else {
            new_seed_stack.push_back(Index_ijk(i-1,j1,k));
          }
        }
        else if (spanw && i>c.si() && (rgnid[i-1][j1][k]!=0 ||
                !(((rgn_inc < 0) && (c[i-1][j1][k] <= c_criterion)) || 
                 ((rgn_inc > 0) && (c[i-1][j1][k] > c_criterion))))) {
          spanw = 0;
        }
        if ((!spane && i<c.ei()) && (rgnid[i+1][j1][k]==0)) {
          if (((rgn_inc < 0) && (c[i+1][j1][k] <= c_criterion)) ||  
              ((rgn_inc > 0) && (c[i+1][j1][k] > c_criterion))) {
            same_fill_stack.push_back(Index_ijk(i+1,j1,k));
            spane = 1;
          }
          else {
            new_seed_stack.push_back(Index_ijk(i+1,j1,k));
            }
        }
        else if (spane && i<c.ei() && (rgnid[i+1][j1][k]!=0 ||
                !(((rgn_inc < 0) && (c[i+1][j1][k] <= c_criterion)) || 
                 ((rgn_inc > 0) && (c[i+1][j1][k] > c_criterion))))) {
          spane = 0;
        }
        if ((!spanb && k>c.sk()) && (rgnid[i][j1][k-1]==0)) {
          if (((rgn_inc < 0) && (c[i][j1][k-1] <= c_criterion)) ||  
              ((rgn_inc > 0) && (c[i][j1][k-1] > c_criterion))) {
            same_fill_stack.push_back(Index_ijk(i,j1,k-1));
            spanb = 1;
          }
          else {
            new_seed_stack.push_back(Index_ijk(i,j1,k-1));
          }
        }
        else if (spanb && k>c.sk() && (rgnid[i][j1][k-1]!=0 ||
                !(((rgn_inc < 0) && (c[i][j1][k-1] <= c_criterion)) || 
                 ((rgn_inc > 0) && (c[i][j1][k-1] > c_criterion))))) {
          spanb = 0;
        }
        if ((!spant && k<c.ek()) && (rgnid[i][j1][k+1]==0)) {
          if (((rgn_inc < 0) && (c[i][j1][k+1] <= c_criterion)) || 
              ((rgn_inc > 0) && (c[i][j1][k+1] > c_criterion))) {
            same_fill_stack.push_back(Index_ijk(i,j1,k+1));
            spant = 1;
          }
          else {
            new_seed_stack.push_back(Index_ijk(i,j1,k+1));
          }
        }
        else if (spant && k<c.ek() && (rgnid[i][j1][k+1]!=0 ||
                  !(((rgn_inc < 0) && (c[i][j1][k+1] <= c_criterion)) || 
                    ((rgn_inc > 0) && (c[i][j1][k+1] > c_criterion))))) {
          spant = 0;
        }
        j1++;
      }          
      else {               
        new_seed_stack.push_back(Index_ijk(i,j1,k));
        break;
      }
    }
  }
}
