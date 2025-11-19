#include "rack.h"

void add_location(const Domain & d,
                  const Range<int> ri,const Range<int> rj,const Range<int> rk,
                  std::vector<Location *> * mons);

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n   - rack's name (it is printed before the values),
*  \param d   - Domain on which the rack is created,
*  \param ri  - range of i's (minimum and maximum)
*  \param j,k - logical coordinates defining racks's position.
*
*  \note This constructor creates the rack in "i" direction.                     
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d, 
           const Range<int> ri, const int j, const int k) 
 : name(n), r_i(ri), r_j(j,j), r_k(k,k) {

  add_location(d,r_i,r_j,r_k,&mons);
#if 0
  int size = ri.last() - ri.first() + 1;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) {
    mons[c++] = new Location(NULL, d, i, j, k);
  }
#endif
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n   - rack's name (it is printed before the values),
*  \param d   - Domain on which the rack is created,
*  \param rj  - range of j's (minimum and maximum)
*  \param i,k - logical coordinates defining racks's position.
*
*  \note This constructor creates the rack in "j" direction.                     
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d, 
           const int i, const Range<int> rj, const int k) 
 : name(n), r_i(i,i), r_j(rj), r_k(k,k) {

  add_location(d,r_i,r_j,r_k,&mons);
#if 0
  int size = rj.last() - rj.first() + 1;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int j=r_j.first(); j<=r_j.last(); j++) {
    mons[c++] = new Location(NULL, d, i, j, k);
  }
#endif
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n   - rack's name (it is printed before the values),
*  \param d   - Domain on which the rack is created,
*  \param rk  - range of k's (minimum and maximum)
*  \param i,j - logical coordinates defining racks's position.
*
*  \note This constructor creates the rack in "k" direction.                     
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d, 
           const int i, const int j, const Range<int> rk) 
 : name(n), r_i(i,i), r_j(j,j), r_k(rk) {

  add_location(d,r_i,r_j,r_k,&mons);
#if 0
  int size = rk.last() - rk.first() + 1;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int k=r_k.first(); k<=r_k.last(); k++) {
    mons[c++] = new Location(NULL, d, i, j, k);
  }
#endif
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n     - rack's name (it is printed before the values),
*  \param d     - Domain on which the rack is created,
*  \param i     - logical coordinates defining racks's position.
*  \param rj,rk - range of j's and k's (minimum and maximum)
*
*  \note This constructor creates the rack in "j-k" plane.
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d,
           const int i, const Range<int> rj, const Range<int> rk)
 : name(n), r_i(i,i), r_j(rj), r_k(rk) {

  //std::cout<<"Rack:begin:"<<boil::cart.iam()<<"\n";
#if 1
  add_location(d,r_i,r_j,r_k,&mons);
#else
  int size_j = rj.last() - rj.first() + 1;
  int size_k = rk.last() - rk.first() + 1;
  int size = size_j*size_k;

  assert(size > 0);

  //std::cout<<"Raci:mons.size()= "<<mons.size()<<"\n";
  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=0; // starts from 1
  for(int j=r_j.first(); j<=r_j.last(); j++) {
    std::cout<<"Rack:j= "<<j<<" iam= "<<boil::cart.iam()<<"\n";
    for(int k=r_k.first(); k<=r_k.last(); k++) {
      //std::cout<<"Rack:j= "<<j<<" k= "<<k<<"\n";
      std::cout<<"c="<<c<<"\n";
      mons[c++] = new Location(NULL, d, i, j, k);
      std::cout<<"c="<<c<<"\n";
      //exit(0);
    }
  }
#endif
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n     - rack's name (it is printed before the values),
*  \param d     - Domain on which the rack is created,
*  \param j     - logical coordinates defining racks's position.
*  \param ri,rk - range of i's and k's (minimum and maximum)
*
*  \note This constructor creates the rack in "i-k" plane.
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d,
           const Range<int> ri, const int j, const Range<int> rk)
 : name(n), r_i(ri), r_j(j,j), r_k(rk) {

  add_location(d,r_i,r_j,r_k,&mons);
#if 0
  int size_i = ri.last() - ri.first() + 1;
  int size_k = rk.last() - rk.first() + 1;
  int size = size_i*size_k;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) {
    for(int k=r_k.first(); k<=r_k.last(); k++) {
      mons[c++] = new Location(NULL, d, i, j, k);
    }
  }
#endif
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param n     - rack's name (it is printed before the values),
*  \param d     - Domain on which the rack is created,
*  \param k     - logical coordinates defining racks's position.
*  \param ri,rj - range of i's and j's (minimum and maximum)
*
*  \note This constructor creates the rack in "i-j" plane.
*******************************************************************************/
Rack::Rack(const char * n, const Domain & d,
           const Range<int> ri, const Range<int> rj, const int k)
 : name(n), r_i(ri), r_j(rj), r_k(k,k) {

  add_location(d,r_i,r_j,r_k,&mons);
#if 0
  int size_i = ri.last() - ri.first() + 1;
  int size_j = rj.last() - rj.first() + 1;
  int size = size_i*size_j;

  assert(size > 0);

  mons.resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) {
    for(int j=r_j.first(); j<=r_j.last(); j++) {
      mons[c++] = new Location(NULL, d, i, j, k);
    }
  }
#endif
}

/******************************************************************************/
void Rack::print(const Scalar & phi) {

  if(name)
    boil::oout << "Rack:" << name << boil::endl;

  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) 
    for(int j=r_j.first(); j<=r_j.last(); j++) 
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        (mons[c++])->print(phi);
      }
}

/******************************************************************************/
void Rack::print(const Vector & u, const Comp & m) {

  if(name)
    boil::oout << "Rack:" << name << boil::endl;

  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++) 
    for(int j=r_j.first(); j<=r_j.last(); j++) 
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        (mons[c++])->print(u,m);
      }
}

/******************************************************************************/
void Rack::save(const Scalar & phi, const char * nm, const int it) {

  std::unique_ptr<real[]> tempArray(new real[mons.size()]); // start from 1

  /* data */
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray[c] = (mons[c])->get_scalar(phi);
        c++;
     }

  if(boil::cart.iam()==0) {

    /* file name */
    std::string fname = name_file(nm, ".mval", it);
    std::cout<<"# Plotting: "<<fname<<"\n";

    /* open a file */
    std::ofstream out(fname.c_str());

    /* header */
    out<<" ist= "<<r_i.first()<<" ied= "<<r_i.last()
       <<" jst= "<<r_j.first()<<" jed= "<<r_j.last()
       <<" kst= "<<r_k.first()<<" ked= "<<r_k.last();
    if(name) out << " Rack= " << name;
    out<<"\n";
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray[c]<<" ";
    }
    out.close();
  }
  //std::cout<<"exit\n";
  //exit(0);

}

/******************************************************************************/
void Rack::save(const Vector & u, const Comp & m,
                const char * nm, const int it) {

  std::unique_ptr<real[]> tempArray(new real[mons.size()]); // start from 1

  /* data */
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray[c] = (mons[c])->get_vector(u,m);
        c++;
     }

  if(boil::cart.iam()==0) {
    /* file name */
    std::string fname = name_file(nm, ".mval", it);

    /* open a file */
    std::ofstream out(fname.c_str());
    std::cout<<"# Plotting: "<<fname<<"\n";

    /* header */
    out<<" ist= "<<r_i.first()<<" ied= "<<r_i.last()
       <<" jst= "<<r_j.first()<<" jed= "<<r_j.last()
       <<" kst= "<<r_k.first()<<" ked= "<<r_k.last();
    if(name) out << " Rack= " << name;
    out<<"\n";

    /* data */
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray[c]<<" ";
    }
    out<<"\n";
    out.close();
  }
}

/******************************************************************************/
void Rack::save_grid(const Scalar & phi, const char * nm) {

  std::unique_ptr<real[]> tempArray_x(new real[mons.size()]); // start from 1
  std::unique_ptr<real[]> tempArray_y(new real[mons.size()]); // start from 1
  std::unique_ptr<real[]> tempArray_z(new real[mons.size()]); // start from 1

  /* data */
  Comp n = Comp::x();
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray_x[c] = (mons[c])->get_grid(phi,n);
        c++;
     }
  n = Comp::y();
  c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray_y[c] = (mons[c])->get_grid(phi,n);
        c++;
     }
  n = Comp::z();
  c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray_z[c] = (mons[c])->get_grid(phi,n);
        c++;
     }

  if(boil::cart.iam()==0) {
    /* file name */
    std::string fname = std::string(nm) +".mgrd";
    //std::string fname = name_file(nm, ".mgrd", 0, boil::cart.iam());
    std::cout<<"# Plotting: "<<fname<<"\n";

    /* open a file */
    std::ofstream out(fname.c_str());

    /* header */
    out<<" ist= "<<r_i.first()<<" ied= "<<r_i.last()
       <<" jst= "<<r_j.first()<<" jed= "<<r_j.last()
       <<" kst= "<<r_k.first()<<" ked= "<<r_k.last();
    if(name) out << " Rack= " << name;
    out<<"\n";

    /* data */
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray_x[c]<<" ";
    }
    out<<"\n";
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray_y[c]<<" ";
    }
    out<<"\n";
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray_z[c]<<" ";
    }
    out<<"\n";
    out.close();
  }
}

/******************************************************************************/
void Rack::save_grid(const Vector & u, const Comp & m, const char * nm) {

  std::unique_ptr<real[]> tempArray_x(new real[mons.size()]); // start from 1
  std::unique_ptr<real[]> tempArray_y(new real[mons.size()]); // start from 1
  std::unique_ptr<real[]> tempArray_z(new real[mons.size()]); // start from 1

  /* data */
  Comp n = Comp::x();
  int c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray_x[c] = (mons[c])->get_grid(u,m,n);
        c++;
     }
  n = Comp::y();
  c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray_y[c] = (mons[c])->get_grid(u,m,n);
        c++;
     }
  n = Comp::z();
  c=1; // starts from 1
  for(int i=r_i.first(); i<=r_i.last(); i++)
    for(int j=r_j.first(); j<=r_j.last(); j++)
      for(int k=r_k.first(); k<=r_k.last(); k++) {
        tempArray_z[c] = (mons[c])->get_grid(u,m,n);
        c++;
     }

  if(boil::cart.iam()==0) {
    /* file name */
    std::string fname = std::string(nm) +".mgrd";
    std::cout<<"# Plotting: "<<fname<<"\n";

    /* open a file */
    std::ofstream out(fname.c_str());

    /* header */
    out<<" ist= "<<r_i.first()<<" ied= "<<r_i.last()
       <<" jst= "<<r_j.first()<<" jed= "<<r_j.last()
       <<" kst= "<<r_k.first()<<" ked= "<<r_k.last();
    if(name) out << " Rack= " << name;
    out<<"\n";

    /* data */
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray_x[c]<<" ";
    }
    out<<"\n";
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray_y[c]<<" ";
    }
    out<<"\n";
    for(int c=1; c<=mons.size()-1; c++) {
      out<<tempArray_z[c]<<" ";
    }
    out<<"\n";
    out.close();
  }
}

/***************************************************************************//**
*  \brief Constructor for monitoring Rack (of points) in computational domain
*
*  \param d         - Domain on which the rack is created,
*  \param ri,rj,rk  - ranges (minimum and maximum)
*******************************************************************************/
void add_location(const Domain & d,
                  const Range<int> ri,const Range<int> rj,const Range<int> rk,
                  std::vector<Location *> * mons){

  int size_i = ri.last() - ri.first() + 1;
  int size_j = rj.last() - rj.first() + 1;
  int size_k = rk.last() - rk.first() + 1;

  int size = size_i * size_j * size_k;

  assert(size > 0);

  mons->resize( size+1 ); // starts from 1

  /* create monitoring points */
  int c=1; // starts from 1
  for(int i=ri.first(); i<=ri.last(); i++)
    for(int j=rj.first(); j<=rj.last(); j++)
      for(int k=rk.first(); k<=rk.last(); k++) {
        (*mons)[c++] = new Location(NULL, d, i, j, k);
      }
}

