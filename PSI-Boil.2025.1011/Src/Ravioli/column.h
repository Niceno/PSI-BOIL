#ifndef COLUMN_H
#define COLUMN_H

/***************************************************************************//**
*  \brief Ravioli class for safe parameter passing to property look-up table.
*******************************************************************************/

//////////////
//          //
//  Column  // 
//          //
//////////////
class Column {
  public:
    explicit Column(const int & v) : val(v)  {};
    operator int () const {return val;}

  private:
    explicit Column() {};
    int val;
};

#endif
