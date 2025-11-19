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
