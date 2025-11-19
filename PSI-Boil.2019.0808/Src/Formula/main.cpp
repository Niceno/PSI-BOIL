#include "formula.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
main()
{
 string request;
 string naredba;

 Formula F;

 cout << "------------------------------------------------------------" << endl;
 cout << " Formula                            Autor : Bojan Niceno    " << endl;
 cout << "------------------------------------------------------------" << endl;
 for(;;)
  {
   cout << "-> ";
   getline(cin, naredba);

//--------------------------------------//
// Naredbe - olitiga rezervirane rijeci //
//--------------------------------------//
   if( naredba == "mem" )
     F.list();
   else if( naredba == "bye" )
     break;
   else
    {
     request = naredba;
     cout << "=> " << F.evaluate(request) << endl;
    }
  }
}
