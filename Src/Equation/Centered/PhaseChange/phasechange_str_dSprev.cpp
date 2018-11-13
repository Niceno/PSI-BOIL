#include "phasechange.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange::str_dSprev() {

#ifdef DEBUG
  std::cout<<"str_dSprev: "<<boil::cart.iam()<<"\n";
#endif

  int kk;
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    nucl->dSprev[i][j] = nucl->area_vapor(i,j,k,Dir::kmin());
    kk=k;
  }
  nucl->store_dSprev=true;

#if 0
  /* debug: output area data to tecplot */
  std::ofstream fout;
  std::string name = "dSprev.dat";
  fout.open(name.c_str());
  fout<<"VARIABLES=\"X\" \"Y\" \"Z\" \"A\" " << boil::endl;
  fout<<"ZONE I="<<dom->ni()-1<<", J="<<dom->nj()-1<<", K="<<2<<"\n";
  fout<<"DATAPACKING=BLOCK, VARLOCATION=([4]=CELLCENTERED)"<<"\n";
  fout << boil::endl;

  for_m(m) {
    if(m==Comp::u()) fout << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) fout << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) fout << "# Z-COORDINATES" << boil::endl;
    int count=0;
    for(int k=kk; k<=kk+1; k++)
      for(int j=1; j<=dom->nj()-1; j++)
        for(int i=1; i<=dom->ni()-1; i++) {
          if(m==Comp::u()) fout << dom->xn(i) << " ";
          if(m==Comp::v()) fout << dom->yn(j) << " ";
          if(m==Comp::w()) fout << dom->zn(k) << " ";
          count++;
          if(count % 8 == 0) fout << boil::endl;
        }
      /* start each coordinate in the new line (if needed) */
      if(count % 8 != 0) fout << boil::endl;
  }

  fout << "# VARIABLE-A" << boil::endl;
  int count=0;
    for(int j=1; j<=dom->nj()-2; j++) {
      for(int i=1; i<=dom->ni()-2; i++) {
        fout << nucl->dSprev[i][j] << " ";
        count++;
        if(count % 8 == 0) fout << boil::endl;
      }
      if(count % 8 != 0) fout << boil::endl;
    }
  fout.close();
  exit(0);
#endif

}

