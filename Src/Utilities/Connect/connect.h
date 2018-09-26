#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <climits>    /* INT_MAX, INT_MIN */

#include "../../Global/global_constants.h"
#include "../../Global/global_precision.h"
#include "../../Global/global_name_file.h"
#include "../../Global/global_malloc.h"
#include "../../Ravioli/range.h"

#ifndef CONNECT_H
#define CONNECT_H

void connect_gmv(const char * bname, const int n, const int t);
void connect_tec(const char * bname, const int n, const int t);
void connect_vtk(const char * bname, const int n, const int t);
void connect_ib_vtk(const char * bname, const int n, const int t);

const real small = boil::nano;

////////////
//        //
//  Node  //
//        //
////////////
class Node {
  public:
    Node(int n_, real x_, real y_, real z_) 
     {x=x_; y=y_; z=z_; old_n=n_; new_n=-1;}

    void new_number(int n_) {new_n = n_;}
    int  new_number() const {return new_n;}
    int  old_number() const {return old_n;}

    friend std::ostream &
      operator << (std::ostream & os, const Node & nod);

    bool operator < (const Node & other) const {
      /*----------+
      |  check x  |
      +----------*/
      if     ( x < other.x-small ) return true;
      else if( x > other.x+small ) return false;
      else { /*----------+
             |  check y  |
             +----------*/
             if     ( y < other.y-small ) return true;
             else if( y > other.y+small ) return false;
             else { /*----------+
                    |  check z  |
                    +----------*/
                    if     ( z < other.z-small ) return true;
                    else                         return false;
        }
      }
    }

  private:
    int  old_n; /* number */
    int  new_n; /* number */
    real x;
    real y;
    real z;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: connect.h,v 1.8 2013/08/23 11:55:26 niceno Exp $'/
+-----------------------------------------------------------------------------*/
