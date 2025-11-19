#ifndef CUT_CELL_H
#define CUT_CELL_H

///////////////
//           //
//  CutCell  //
//           //
///////////////
class CutCell {

  public:

    //! For setting the geometrical values.
    void fS(const real val, const int d, const int o) {f_s[d][o]=val;}
    void fV(const real val)                           {f_v=val;}
    void fdxw(const real val)                         {f_dxw=val;}
    void fdxe(const real val)                         {f_dxe=val;}
    void fdys(const real val)                         {f_dys=val;}
    void fdyn(const real val)                         {f_dyn=val;}
    void fdzb(const real val)                         {f_dzb=val;}
    void fdzt(const real val)                         {f_dzt=val;}
    void fE(const real val, const int d, const int i1, const int i2)
                                                      {f_e[d][i1][i2]=val;}
    void fP(const real val, const int ii, const int jj, const int kk)
                                                      {f_p[ii][jj][kk]=val;}
    void set_iacpt(const int val)                     {i_acpt=val;}
    void set_jacpt(const int val)                     {j_acpt=val;}
    void set_kacpt(const int val)                     {k_acpt=val;}

    //! For using the geometrical values.
    real fS(const int d, const int o) const {return f_s[d][o];}
    real fV()                         const {return f_v;}
    real fdxw()                       const {return f_dxw;}
    real fdxe()                       const {return f_dxe;}
    real fdys()                       const {return f_dys;}
    real fdyn()                       const {return f_dyn;}
    real fdzb()                       const {return f_dzb;}
    real fdzt()                       const {return f_dzt;}
    real fE(const int d, const int i1, const int i2)
                                      const {return f_e[d][i1][i2];}
    real fP(const int ii, const int jj, const int kk) 
                                      const {return f_p[ii][jj][kk];}
    int iacpt()                       const {return i_acpt;}
    int jacpt()                       const {return j_acpt;}
    int kacpt()                       const {return k_acpt;}

    //! Set the position of the cell.
    void ijk(const int i_, const int j_, const int k_) {
      i = i_;  j = j_;  k = k_;
    }

    //! Get the position of the cell.
    void ijk(int * i_, int * j_, int * k_) const {
      *i_ = i;  *j_ = j;  *k_ = k;
    }

    //int nccells_in,nccells_bd; // no. cell of inside and boundary

  private:
    real f_s[3][2];
    real f_v;
    real f_dxw, f_dxe, f_dys, f_dyn, f_dzb, f_dzt;
    real f_e[3][2][2];
    real f_p[2][2][2];
    int  i,j,k;
    int  i_acpt, j_acpt, k_acpt;
};

#endif
