#include <iostream>
#include <cassert>
#include <cmath>  
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#define real double

///////////
//       //
//  Get  //
//       //
///////////
class Get {
  public:
    static const Get rho   () {return Get(2);}
    static const Get mu    () {return Get(3);}
    static const Get cp    () {return Get(4);}
    static const Get lambda() {return Get(5);}
    static const Get gamma () {return Get(6);}
    static const Get sigma () {return Get(7);}

    operator int () const {return val;}

  private:
    /* prevent creation of new options */
    explicit Get(const bool i) {val = i;}
    int val;
};

/*============================================================================*/
real get_property(const std::vector < std::vector < real > > t, 
                  const real val,
                  const Get & g) {
/*

   t:   0     1     2     3     4     5     6     7
        |-----|-----|-----|-----|-----|-----|-----|
   seg:    1     2     3     4     5     6     7

*/

  assert(val >= t[0][0]);
  assert(val <= t[t.size()-1][0]);

  int lowseg  = 1;
  int highseg = t.size()-1;

  int currseg = (highseg+lowseg)/2;

  for(;;) {
    if(val > t[currseg][0]) {
      lowseg=currseg;
      currseg = (int) ceil(((real)highseg+(real)lowseg)/2.0);
    } else if(val < t[currseg-1][0]) {
      highseg=currseg;
      currseg = (int) floor(((real)highseg+(real)lowseg)/2.0);
    } else {
      std::cout << "low  = " << t[currseg-1][0] << std::endl;
      std::cout << "val  = " << val             << std::endl;
      std::cout << "high = " << t[currseg]  [0] << std::endl;
      const real high = t[currseg]  [0];
      const real  low = t[currseg-1][0];
      const real hcoef = (val -low) / (high-low);
      const real lcoef = (high-val) / (high-low);
      if( g==Get::rho() ) 
        return lcoef * t[currseg-1][3] + hcoef * t[currseg][3];
      else if( g==Get::rho() ) 
        return lcoef * t[currseg-1][3] + hcoef * t[currseg][3];
      else if( g==Get::cp() ) 
        return lcoef * t[currseg-1][3] + hcoef * t[currseg][3];
      else if( g==Get::lambda() ) 
        return lcoef * t[currseg-1][3] + hcoef * t[currseg][3];
      else if( g==Get::mu() ) 
        return lcoef * t[currseg-1][3] + hcoef * t[currseg][3];
    }
  }

  return -1;
}

/*============================================================================*/
int main() {

  std::ifstream in;
  in.open("co2.txt");

  /* stop if file is not present */
  if( in.rdstate() != 0 ) {
    std::cout << "failed to open co2.txt; " << "exiting!" << std::endl;
    exit(0);
  }

  std::vector < std::vector < real > > table;
  std::string line;

  /* skip first two lines */
  getline(in, line);
  getline(in, line);

  /* read the rest of the lines */
  while ( getline(in, line) ) {
    std::vector < real > data;
    real value;
    std::istringstream iss(line);
    while (iss >> value)
     {
      data.push_back(value);
     }
    table.push_back(data);
  }

  std::cout << "table size: " << table.size() << std::endl;

  std::cout << "property = " << get_property(table, 57.59999, Get::rho()) << std::endl;

}
