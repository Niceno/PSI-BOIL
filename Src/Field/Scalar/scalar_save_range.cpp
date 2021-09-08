#include "scalar.h"
//#include "../../Domain/domain.h"

// Global I,J,K are used for the intex
/******************************************************************************/
void Scalar::save_range(const Range<int> RI, const int J, const int K,
                        const char * nam, const int it) {
  const Range<int> RJ(J,J);
  const Range<int> RK(K,K);
  save_range(RI,RJ,RK,nam,it);
}
/******************************************************************************/
void Scalar::save_range(const int I, const Range<int> RJ, const int K,
                        const char * nam, const int it) {
  const Range<int> RI(I,I);
  const Range<int> RK(K,K);
  save_range(RI,RJ,RK,nam,it);
}
/******************************************************************************/
void Scalar::save_range(const int I, const int J, const Range<int> RK,
                        const char * nam, const int it) {
  const Range<int> RI(I,I);
  const Range<int> RJ(J,J);
  save_range(RI,RJ,RK,nam,it);
}

/******************************************************************************/
void Scalar::save_range(const int I, const Range<int> RJ, const Range<int> RK,
                        const char * nam, const int it) {
  const Range<int> RI(I,I);
  save_range(RI,RJ,RK,nam,it);
}
/******************************************************************************/
void Scalar::save_range(const Range<int> RI, const int J, const Range<int> RK,
                        const char * nam, const int it) {
  const Range<int> RJ(J,J);
  save_range(RI,RJ,RK,nam,it);
}
/******************************************************************************/
void Scalar::save_range(const Range<int> RI, const Range<int> RJ, const int K,
                        const char * nam, const int it) {
  const Range<int> RK(K,K);
  save_range(RI,RJ,RK,nam,it);
}

/******************************************************************************/
void Scalar::save_range(const Range<int> RI, const Range<int> RJ,
                        const Range<int> RK,
                        const char * nam, const int it) {
  /* file name */
  std::string name = name_file(nam, ".dat", it);
  boil::oout<<"# Plotting: "<<name<<"\n";

  /* open a file */
  std::ofstream out(name.c_str());
  if(boil::cart.iam()==0) {
    out<<"VARIABLES = \"X\" \"Y\" \"Z\" \"A\"\n";
    out<<"ZONE, I="<<RI.last()-RI.first()+1
       <<", J="<<RJ.last()-RJ.first()+1
       <<", K="<<RK.last()-RK.first()+1<<", DATAPACKING=POINT\n";
  }

  //for(int I=RI.first(); I<=RI.last(); I++) {
  //  for(int J=RJ.first(); J<=RJ.last(); J++) {
  //    for(int K=RK.first(); K<=RK.last(); K++) {
  for(int K=RK.first(); K<=RK.last(); K++) {
    for(int J=RJ.first(); J<=RJ.last(); J++) {
      for(int I=RI.first(); I<=RI.last(); I++) {
        real x=0.0, y=0.0, z=0.0;
        real v=0.0;
        if( dom->contains_IJK(I,J,K) ) {
          int i=I, j=J, k=K;
          dom->locals(&i,&j,&k); // change i,j,k
          x=xc(i);
          y=yc(j);
          z=zc(k);
          v=val[i][j][k];   // with changed i,j,k
        }
        //real v = (*this)(I,J,K);
        boil::cart.sum_real(&x);
        boil::cart.sum_real(&y);
        boil::cart.sum_real(&z);
        boil::cart.sum_real(&v);
        if(boil::cart.iam()==0)
          out<<x<<" "<<y<<" "<<z<<" "<<v<<"\n";
      }
    }
  }

  /* close a file */
  out.close();
}
