#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::plot_f(const char * nam) {
/***************************************************************************//**
*  \brief Output for TECPLOT on node.
*******************************************************************************/
  std::ofstream fout;
  std::string name = name_file(nam, ".dat", 3, boil::cart.iam());
  fout.open(name.c_str());
  fout<<"VARIABLES=\"X\" \"Y\" \"Z\" \"A\" " << boil::endl;
  fout<<"ZONE I="<<ei()-si()+2<<", J="
                 <<ej()-sj()+2<<", K="
                 <<ek()-sk()+2<<"\n";
  fout<<"DATAPACKING=BLOCK, VARLOCATION=([4]=NODAL)"<<"\n";
  fout << boil::endl;
  /* EXTENDED BUFFERS HINT:
     Line 12 should probably be replaced with this:
  fout<<"ZONE I="<<dom->ei()-dom->si()+1 <<", 
              J="<<dom->ej()-dom->sj()+1 <<", 
              K="<<dom->ek()-dom->sk()+1<<"\n";
  */

  for_m(m) {
    if(m==Comp::u()) fout << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) fout << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) fout << "# Z-COORDINATES" << boil::endl;
    int count=0;
    for(int k=sk(); k<=ek()+1; k++)
      for(int j=sj(); j<=ej()+1; j++)
        for(int i=si(); i<=ei()+1; i++) {
          if(m==Comp::u()) fout << dom->xn(i) << " ";
          if(m==Comp::v()) fout << dom->yn(j) << " ";
          if(m==Comp::w()) fout << dom->zn(k) << " ";
          count++;
          if(count % 8 == 0) fout << boil::endl;
        }
     /* start each coordinate in the new line (if needed) */
     if(count % 8 != 0) fout << boil::endl;
   }

  int count=0;
  for(int k=sk(); k<=ek()+1; k++)
    for(int j=sj(); j<=ej()+1; j++)
      for(int i=si(); i<=ei()+1; i++) {
        fout << scheme.f[i][j][k] << " ";
        count++;
        if(count % 8 == 0) fout << boil::endl;
      }
    /* EXTENDED BUFFERS HINT:
       Lines 47-49 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek()+1; k++)
      for(int j=sj(); j<=dom->ej()+1; j++)
        for(int i=si(); i<=dom->ei()+1; i++) {
    */
}

/******************************************************************************/
void CIPCSL2::plot_sigx(const char * nam) {
/***************************************************************************//**
*  \brief Output for TECPLOT on edge-i.
*******************************************************************************/
  std::ofstream fout;
  std::string name = name_file(nam, ".dat", 3, boil::cart.iam());
  fout.open(name.c_str());
  fout<<"VARIABLES=\"X\" \"Y\" \"Z\" \"A\" " << boil::endl;
  fout<<"ZONE I="<<ei()-si()+1
      <<", J="<<ej()-sj()+2
      <<", K="<<ek()-sk()+2<<"\n";
  fout<<"DATAPACKING=BLOCK, VARLOCATION=([4]=NODAL)"<<"\n";
  fout << boil::endl;
  /* EXTENDED BUFFERS HINT:
     Line 71 should probably be replaced with this:
  fout<<"ZONE I="<<dom->ei()-dom->si()   <<", 
              J="<<dom->ej()-dom->sj()+1 <<", 
              K="<<dom->ek()-dom->sk()+1<<"\n";
  */

  for_m(m) {
    if(m==Comp::u()) fout << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) fout << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) fout << "# Z-COORDINATES" << boil::endl;
    int count=0;
    for(int k=sk(); k<=ek()+1; k++)
      for(int j=sj(); j<=ej()+1; j++)
        for(int i=si(); i<=ei(); i++) {
          if(m==Comp::u()) fout << dom->xc(i) << " ";
          if(m==Comp::v()) fout << dom->yn(j) << " ";
          if(m==Comp::w()) fout << dom->zn(k) << " ";
          count++;
          if(count % 8 == 0) fout << boil::endl;
        }
     /* start each coordinate in the new line (if needed) */
     if(count % 8 != 0) fout << boil::endl;
   }
    /* EXTENDED BUFFERS HINT:
       Lines 86-88 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek()+1; k++)
      for(int j=sj(); j<=dom->ej()+1; j++)
        for(int i=si(); i<=dom->ei(); i++) {
    */

  int count=0;
  for(int k=sk(); k<=ek()+1; k++)
    for(int j=sj(); j<=ej()+1; j++)
      for(int i=si(); i<=ei(); i++) {
        real dx=phi.dxc(i);
        fout << scheme.sigx[i][j][k]/dx << " ";
        count++;
        if(count % 8 == 0) fout << boil::endl;
      }
    /* EXTENDED BUFFERS HINT:
       Lines 106-108 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek()+1; k++)
      for(int j=sj(); j<=dom->ej()+1; j++)
        for(int i=si(); i<=dom->ei(); i++) {
    */
}


/******************************************************************************/
void CIPCSL2::plot_sigy(const char * nam) {
/***************************************************************************//**
*  \brief Output for TECPLOT on edge-j.
*******************************************************************************/
  std::ofstream fout;
  std::string name = name_file(nam, ".dat", 3, boil::cart.iam());
  fout.open(name.c_str());
  fout<<"VARIABLES=\"X\" \"Y\" \"Z\" \"A\" " << boil::endl;
  fout<<"ZONE I="<<ei()-si()+2<<", J="
                 <<ej()-sj()+1<<", K="
                 <<ek()-sk()+2<<"\n";
  fout<<"DATAPACKING=BLOCK, VARLOCATION=([4]=NODAL)"<<"\n";
  fout << boil::endl;
  /* EXTENDED BUFFERS HINT:
     Line 132 should probably be replaced with this:
  fout<<"ZONE I="<<dom->ei()-dom->si()+1 <<", 
              J="<<dom->ej()-dom->sj()   <<", 
              K="<<dom->ek()-dom->sk()+1<<"\n";
  */

  for_m(m) {
    if(m==Comp::u()) fout << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) fout << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) fout << "# Z-COORDINATES" << boil::endl;
    int count=0;
    for(int k=sk(); k<=ek()+1; k++)
      for(int j=sj(); j<=ej()  ; j++)
        for(int i=si(); i<=ei()+1; i++) {
          if(m==Comp::u()) fout << dom->xn(i) << " ";
          if(m==Comp::v()) fout << dom->yc(j) << " ";
          if(m==Comp::w()) fout << dom->zn(k) << " ";
          count++;
          if(count % 8 == 0) fout << boil::endl;
        }
     /* start each coordinate in the new line (if needed) */
     if(count % 8 != 0) fout << boil::endl;
   }
    /* EXTENDED BUFFERS HINT:
       Lines 150-152 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek()+1; k++)
      for(int j=sj(); j<=dom->ej()  ; j++)
        for(int i=si(); i<=dom->ei()+1; i++) {
    */

  int count=0;
  for(int k=sk(); k<=ek()+1; k++)
    for(int j=sj(); j<=ej()  ; j++)
      for(int i=si(); i<=ei()+1; i++) {
        real dy=phi.dyc(j);
        fout << scheme.sigy[i][j][k]/dy << " ";
        count++;
        if(count % 8 == 0) fout << boil::endl;
      }
    /* EXTENDED BUFFERS HINT:
       Lines 167-169 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek()+1; k++)
      for(int j=sj(); j<=dom->ej()  ; j++)
        for(int i=si(); i<=dom->ei()+1; i++) {
    */
}

/******************************************************************************/
void CIPCSL2::plot_sigz(const char * nam) {
/***************************************************************************//**
*  \brief Output for TECPLOT on edge-k.
*******************************************************************************/
  std::ofstream fout;
  std::string name = name_file(nam, ".dat", 3, boil::cart.iam());
  fout.open(name.c_str());
  fout<<"VARIABLES=\"X\" \"Y\" \"Z\" \"A\" " << boil::endl;
  fout<<"ZONE I="<<ei()-si()+2<<", J="
                 <<ej()-sj()+2<<", K="
                 <<ek()-sk()+1<<"\n";
  fout<<"DATAPACKING=BLOCK, VARLOCATION=([4]=NODAL)"<<"\n";
  fout << boil::endl;
  /* EXTENDED BUFFERS HINT:
     Line 192 should probably be replaced with this:
  fout<<"ZONE I="<<dom->ei()-dom->si()+1 <<", 
              J="<<dom->ej()-dom->sj()+1 <<", 
              K="<<dom->ek()-dom->sk()  <<"\n";
  */

  for_m(m) {
    if(m==Comp::u()) fout << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) fout << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) fout << "# Z-COORDINATES" << boil::endl;
    int count=0;
    for(int k=sk(); k<=ek();   k++)
      for(int j=sj(); j<=ej()+1; j++)
        for(int i=si(); i<=ei()+1; i++) {
          if(m==Comp::u()) fout << dom->xn(i) << " ";
          if(m==Comp::v()) fout << dom->yn(j) << " ";
          if(m==Comp::w()) fout << dom->zc(k) << " ";
          count++;
          if(count % 8 == 0) fout << boil::endl;
        }
     /* start each coordinate in the new line (if needed) */
     if(count % 8 != 0) fout << boil::endl;
   }
    /* EXTENDED BUFFERS HINT:
       Lines 207-209 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek();   k++)
      for(int j=sj(); j<=dom->ej()+1; j++)
        for(int i=si(); i<=dom->ei()+1; i++) {
    */

  int count=0;
  for(int k=sk(); k<=ek();   k++)
    for(int j=sj(); j<=ej()+1; j++)
      for(int i=si(); i<=ei()+1; i++) {
        real dz=phi.dzc(k);
        fout << scheme.sigz[i][j][k]/dz << " ";
        count++;
        if(count % 8 == 0) fout << boil::endl;
      }
    /* EXTENDED BUFFERS HINT:
       Lines 227-229 should probably be replaced with this:
    for(int k=sk(); k<=dom->ek();   k++)
      for(int j=sj(); j<=dom->ej()+1; j++)
        for(int i=si(); i<=dom->ei()+1; i++) {
    */
}

/******************************************************************************/
void CIPCSL2::plot_sxyz(const char * nam, const Comp & mc) {
/***************************************************************************//**
*  \brief Output for TECPLOT on face.
*******************************************************************************/
  int ir,jr,kr;
  if(mc==Comp::i()){
    ir=1; jr=0;   kr=0;
  } else if(mc==Comp::j()){
    ir=0;   jr=1; kr=0;
  } else if(mc==Comp::k()){
    ir=0;   jr=0; kr=1;
  }
  /* EXTENDED BUFFERS HINT:
     Lines 249-255 should probably be replaced with this:
  if(mc==Comp::i()){
    ir=2*boil::BW-1; jr=2*boil::BW;   kr=2*boil::BW;
  } else if(mc==Comp::j()){
    ir=2*boil::BW;   jr=2*boil::BW-1; kr=2*boil::BW;
  } else if(mc==Comp::k()){
    ir=2*boil::BW;   jr=2*boil::BW;   kr=2*boil::BW-1;
  */

  std::ofstream fout;
  std::string name = name_file(nam, ".dat", 3, boil::cart.iam());
  fout.open(name.c_str());
  fout<<"VARIABLES=\"X\" \"Y\" \"Z\" \"A\" " << boil::endl;
  fout<<"ZONE I="<<ei()-si()+1+ir
         <<", J="<<ej()-sj()+1+jr
         <<", K="<<ek()-sk()+1+kr<<boil::endl;
  fout<<"DATAPACKING=BLOCK, VARLOCATION=([4]=NODAL)"<<boil::endl;
  fout << boil::endl;

  for_m(m) {
    if(m==Comp::u()) fout << "# X-COORDINATES" << boil::endl;
    if(m==Comp::v()) fout << "# Y-COORDINATES" << boil::endl;
    if(m==Comp::w()) fout << "# Z-COORDINATES" << boil::endl;
    int count=0;
    for(int k=sk(); k<=ek()+kr; k++)
      for(int j=sj(); j<=ej()+jr; j++)
        for(int i=si(); i<=ei()+ir; i++) {
       if(m==Comp::u() && mc==Comp::i()) fout << dom->xn(i) << " ";
       if(m==Comp::u() && mc!=Comp::i()) fout << dom->xc(i) << " ";
       if(m==Comp::v() && mc==Comp::j()) fout << dom->yn(j) << " ";
       if(m==Comp::v() && mc!=Comp::j()) fout << dom->yc(j) << " ";
       if(m==Comp::w() && mc==Comp::k()) fout << dom->zn(k) << " ";
       if(m==Comp::w() && mc!=Comp::k()) fout << dom->zc(k) << " ";
       count++;
       if(count % 8 == 0) fout << boil::endl;
    }
     /* start each coordinate in the new line (if needed) */
     if(count % 8 != 0) fout << boil::endl;
   }
    /* EXTENDED BUFFERS HINT:
       Lines 281-283 should probably be replaced with this:
    for_mijk(mc,i,j,k)
    */

  int count=0;
  for(int k=sk(); k<=ek()+kr; k++)
    for(int j=sj(); j<=ej()+jr; j++)
      for(int i=si(); i<=ei()+ir; i++) {
        real dx=phi.dxc(i);
        real dy=phi.dyc(j);
        real dz=phi.dzc(k);
        real atmp;
        if(mc==Comp::i())atmp=dy*dz;
        if(mc==Comp::j())atmp=dx*dz;
        if(mc==Comp::k())atmp=dx*dy;
        fout << sxyz[mc][i][j][k]/atmp << " ";
        count++;
        if(count % 8 == 0) fout << boil::endl;
  }
}
