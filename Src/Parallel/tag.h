/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing.
*
*  Used inside Communicator to denote the message tag.
*******************************************************************************/
#ifndef TAG_H
#define TAG_H

///////////
//       //
//  Tag  //
//       //
///////////
class Tag {
  public:
    explicit Tag(const int i) : val(i) {}

    operator int() const {return val;} 
  private:
    const int val;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: tag.h,v 1.1 2009/07/17 16:37:18 niceno Exp $'/
+-----------------------------------------------------------------------------*/
