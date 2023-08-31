#include "domain.h"

/******************************************************************************/
void Domain::save(const char * nm) {

  if(boil::cart.iam()==0) {
    /* file name */
    std::stringstream stream;
    stream << nm;
    std::string fname;
    stream >> fname;
    std::cout<<fname<<"\n";

    /* open a file */
    std::ofstream ofs(fname);

    /* write data */
    int n_x = gi();
    int n_y = gj();
    int n_z = gk();
    ofs.write(reinterpret_cast<const char *> (&n_x), sizeof(int));
    ofs.write(reinterpret_cast<const char *> (&n_y), sizeof(int));
    ofs.write(reinterpret_cast<const char *> (&n_z), sizeof(int));
    std::cout<<"Domain::save::ni,nj,nk= "<<n_x<<" "<<n_y<<" "<<n_z<<"\n";

    real tmp;
    for(int i=0; i<n_x; i++){
      tmp = xn_global(i);
      //std::cout<<"i,xn= "<<i<<" "<<tmp<<"\n";
      ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    }
    for(int j=0; j<n_y; j++){
      tmp = yn_global(j);
      //std::cout<<"j,yn= "<<j<<" "<<tmp<<"\n";
      ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    }
    for(int k=0; k<n_z; k++){
      tmp = zn_global(k);
      //std::cout<<"k,zn= "<<k<<" "<<tmp<<"\n";
      ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    }

    //for(int i=0; i<n_x; i++){
    tmp = xc_global(0);
    //std::cout<<"i,xn_stg= "<<0<<" "<<tmp<<"\n";
    ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    for(int i=1; i<n_x+1; i++){
      tmp = xc_global(i-1);
      //std::cout<<"i,xn_stg= "<<i<<" "<<tmp<<"\n";
      ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    }

    tmp = yc_global(0);
    //std::cout<<"j,yn_stg= "<<0<<" "<<tmp<<"\n";
    ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    for(int j=1; j<n_y+1; j++){
      tmp = yc_global(j-1);
      //std::cout<<"j,yn_stg= "<<j<<" "<<tmp<<"\n";
      ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    }

    tmp = zc_global(0);
    //std::cout<<"k,zn_stg= "<<0<<" "<<tmp<<"\n";
    ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    for(int k=1; k<n_z+1; k++){
      tmp = zc_global(k-1);
      //std::cout<<"k,zn_stg= "<<k<<" "<<tmp<<"\n";
      ofs.write(reinterpret_cast<const char *> (&tmp),sizeof(real));
    }

    /* close a file */
    ofs.close();
  }

}
#if 0
/******************************************************************************/
int Domain::I(real x) const {

  for(int i_=0; i_<ni(); i_++) 
    if(x >= xn(i_) && x <= xn(i_+1)) return i_ + cxg().first() - 1;
  return -1;
}

/******************************************************************************/
int Domain::J(real y) const {

  for(int j_=0; j_<nj(); j_++) 
    if(y >= yn(j_) && y <= yn(j_+1)) return j_ + cyg().first() - 1;
  return -1;
}

/******************************************************************************/
int Domain::K(real z) const {

  for(int k_=0; k_<nk(); k_++) 
    if(z >= zn(k_) && z <= zn(k_+1)) return k_ + czg().first() - 1;
  return -1;
}
#endif
