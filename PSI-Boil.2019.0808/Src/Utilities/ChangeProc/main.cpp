#include "changeproc.h"
using namespace std;

/* important: does not work with custom domain decomposition!! */
/* important: factor and decomp must be consistent with the distribute
              member function of Domain class!! */

void factor(int n, int * factor, int * number);
void decomp(int * id, const int i1, const int i2, const int i3, const int i4 );
void scalar( const char *nm, const int nic, const int njc, const int nkc
             , const int ts, const int np, const int np_new
             , const int ii, const int jj, const int kk);
void vector_treat( const char *nm, const int nic, const int njc, const int nkc
             , const int ts, const int np, const int np_new);
void nucleation( const char *nm, const int nic, const int njc, const int nkc
             , const int ts, const int np, const int np_new);

/******************************************************************************/
int main(int argc, char * argv[]) {

  cout<<"#### Change number of processors. ####\n\n";

  cout<<"### Reading input deck. ####\n";

  ifstream input_deck("input.txt");
  istream * inp = &cin;

  if(input_deck.is_open()) {
    inp  = &input_deck;
  } else {
    cout<<"# Input file *input.txt* not found. "
             <<"Defaulting to standard input.\n";
  }

  int nic, njc, nkc, ts, np, np_new;

  cout<<"# Enter nicell, njcell, nkcell.\n";
  *inp>>nic>>njc>>nkc;
  cout<<"# nicell= "<<nic<<" njcell= "<<njc<<" nkcell= "<<nkc<<"\n";

  cout<<"# Enter current np (No. processors).\n";
  *inp>>np;
  cout<<"# np= "<<np<<"\n";

  cout<<"# Enter new np (No. processors).\n";
  *inp>>np_new;
  cout<<"# np_new= "<<np_new<<"\n";

  int ts_s, ts_e, ts_i;
  cout<<"# Enter start, end & interval of time step.\n";
  *inp>>ts_s>>ts_e>>ts_i;
  cout<<"# time step= "<<ts_s<<" "<<ts_e<<" "<<ts_i<<"\n";

  while (true) {

    cout<<"# Enter file name. (Example: press)\n";
    cout<<"# In case of CIPCSL2, type object name (without -f, -sigx, etc.)\n";
    string ss;
    *inp>>ss;
    const char *nm = ss.c_str();
    cout<<"# filename= "<<ss<<"\n";

    /* file type */
    int itype;
    cout<<"# Enter file type.  0:Scalar, 1:Vector, 2:CIPCSL2, 3:Nucleation\n";
    *inp>>itype;
    cout<<"# file type = "<<itype<<"\n";

    cout<<"### Input deck built successfully. ###\n";

    for (int ts = ts_s; ts <= ts_e; ts +=ts_i) {

      cout<<"\n## processing time step "<<ts<<endl;

      string str_tmp;

      switch(itype) {

        case 0 :
          scalar( nm, nic, njc, nkc, ts, np, np_new,0,0,0 );
          break;
        case 1 :
          vector_treat( nm, nic, njc, nkc, ts, np, np_new );
          break;
        case 2 :
          /* clr: cell center */
          str_tmp = nm;
          str_tmp += "-phi";
          scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,0,0,0 );
          /* clr: cell center */
          str_tmp = nm;
          str_tmp += "-clr";
          scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,0,0,0 );
          /* f: node (vertex) */
          str_tmp = nm;
          str_tmp += "-f";
          scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,1,1,1 );
          /* sigx: edge-x */
          str_tmp = nm;
          str_tmp += "-sigx";
          scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,0,1,1 );
          /* sigx: edge-y */
          str_tmp = nm;
          str_tmp += "-sigy";
          scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,1,0,1 );
          /* sigx: edge-y */
          str_tmp = nm;
          str_tmp += "-sigz";
          scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,1,1,0 );
          /* sxyz: vecto */
          str_tmp = nm;
          str_tmp += "-face";
          vector_treat ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new );
          break;
        case 3 :
          nucleation ( nm, nic, njc, nkc, ts, np, np_new );
          break;
        default :
          cout<<" Unrecognized file type. Exiting."<<endl;
          exit(0);
      }
    }
    int ians;
    cout<<"# continue? 1(yes) or 0(no).\n";
    cin>>ians;
    if (ians!=1) break;

  }

  return 0;

}

/*****************************************************************************/
void scalar( const char *nm, const int nic, const int njc, const int nkc, 
             const int ts, const int np, const int np_new, 
             const int ii, const int jj, const int kk){
    /* allocate */
    int ni = nic+2*boil::BW;
    int nj = njc+2*boil::BW;
    int nk = nkc+2*boil::BW;
    real *** q;
    alloc3d( &q, ni+ii, nj+jj, nk+kk );

    for(int i=0; i<ni; i++)
      for(int j=0; j<nj; j++)
        for(int k=0; k<nk; k++)
          q[i][j][k]=1.0e+300;

    string name;
    /* get decomposed size */
    name = name_file(nm, ".bck", ts, 0); 
    ifstream in(name.c_str(), ios::binary);

    /* stop if file is not present */
    if( in.rdstate() != 0 ) {
      cout << "failed to open " << name << endl;
      cout << "exiting!" << endl;
      exit(0);
    }

    int ni_dec, nj_dec, nk_dec;
    in.read(reinterpret_cast<char *> (&ni_dec), sizeof(int));
    in.read(reinterpret_cast<char *> (&nj_dec), sizeof(int));
    in.read(reinterpret_cast<char *> (&nk_dec), sizeof(int));
    cout<<"# ni_dec, nj_dec, nk_dec= "
             <<ni_dec<<" "<<nj_dec<<" "<<nk_dec<<"\n";

    /* allocate decomposed array */
    real *** val;
    alloc3d( &val, ni_dec, nj_dec, nk_dec );

    int nic_dec = ni_dec-2*boil::BW-ii;
    int njc_dec = nj_dec-2*boil::BW-jj;
    int nkc_dec = nk_dec-2*boil::BW-kk;
    cout<<"# nic_dec, njc_dec, nkc_dec= "
             <<nic_dec<<" "<<njc_dec<<" "<<nkc_dec<<"\n";

    int ni_block = nic/nic_dec;
    int nj_block = njc/njc_dec;
    int nk_block = nkc/nkc_dec;
    cout<<"# ni_block, nj_block, nk_block= "
             <<ni_block<<" "<<nj_block<<" "<<nk_block<<"\n";

    /*----------+
    |  compose  |
    +----------*/
    for (int iproc = 0; iproc < np; iproc++) {

      int i_block = real(iproc + 0)/real(nj_block*nk_block) + 1;
      int j_block = real(iproc + 0 - (i_block-1)*(nj_block*nk_block))
                  / real(nk_block) + 1;
      int k_block = real(iproc + 1 - (i_block-1)*(nj_block*nk_block)
                        -(j_block-1)*nk_block); 
      cout<<"iproc= "<<iproc<<" "<<i_block<<" "<<j_block<<" "<<k_block<<"\n";

      /* file name */
      name = name_file(nm, ".bck", ts, iproc);

      /* open a file */
      ifstream in(name.c_str(), ios::binary);

      /* stop if file is not present */
      if( in.rdstate() != 0 ) {
        cout << "failed to open " << name << endl;
        cout << "exiting!" << endl;
        exit(0);
      }

      /* read val */
      int ni_tmp, nj_tmp, nk_tmp;
      in.read(reinterpret_cast<char *> (&ni_tmp), sizeof(int));
      in.read(reinterpret_cast<char *> (&nj_tmp), sizeof(int));
      in.read(reinterpret_cast<char *> (&nk_tmp), sizeof(int));
      in.read(reinterpret_cast<char *> (val[0][0]),
                              ni_dec*nj_dec*nk_dec*sizeof(real));

      /* copy val to q */
      int ist = (i_block-1) * nic_dec;
      int jst = (j_block-1) * njc_dec;
      int kst = (k_block-1) * nkc_dec;
      //cout<<"ist,jst,kst = "<<ist<<" "<<jst<<" "<<kst<<"\n";
      int imax=0; 
      int jmax=0; 
      int kmax=0; 
      int istl=boil::BW; int jstl=boil::BW; int kstl=boil::BW;
      if (i_block==1) istl=0;
      if (j_block==1) jstl=0;
      if (k_block==1) kstl=0;
      for(int i=istl; i<ni_dec; i++)
        for(int j=jstl; j<nj_dec; j++)
          for(int k=kstl; k<nk_dec; k++) {
            q[ist+i][jst+j][kst+k] = val[i][j][k];
            if(ist+i>imax) imax=ist+i;
            if(jst+j>jmax) jmax=jst+j;
            if(kst+k>kmax) kmax=kst+k;
          }
      //cout<<"imax, jmax, kmax= "<<imax<<" "<<jmax<<" "<<kmax<<"\n";

    }

    /*-------------+
    |  deallocate  |
    +-------------*/
    dealloc3d( &val);

    for(int i=0; i<ni; i++)
      for(int j=0; j<nj; j++)
        for(int k=0; k<nk; k++) {
          if(q[i][j][k]>1.0e+299) {
            cout<<"Error! q is not imported.\n";
            exit(0);
          }
        }

    /*----------------------------+
    |  output bck to single-file  |
    +----------------------------*/
    if (np_new==1) {
      /* file name */
      name = name_file(nm, ".bck2", ts, 0);

      int ni_out=ni+ii;
      int nj_out=nj+jj;
      int nk_out=nk+kk;

      /* open a file */
      ofstream out(name.c_str(), ios::binary);

      out.write(reinterpret_cast<const char *> (&ni_out), sizeof(int));
      out.write(reinterpret_cast<const char *> (&nj_out), sizeof(int));
      out.write(reinterpret_cast<const char *> (&nk_out), sizeof(int));
      out.write(reinterpret_cast<const char *> (q[0][0]),
              ni_out*nj_out*nk_out*sizeof(real));
    } else {

      /*------------+
      |  decompose  |
      +------------*/
      cout<<"Gotcha!!!!";
      int * dims = new int[3];
      decomp(dims, np_new, nic, njc, nkc);
      cout<<"# dims[]= "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"\n";
      int nic_dec_new=nic/dims[0];
      int njc_dec_new=njc/dims[1];
      int nkc_dec_new=nkc/dims[2];
      cout<<"# nic_dec_new= "<<nic_dec_new<<" "<<njc_dec_new<<" "
          <<nkc_dec_new<<"\n";
      int ni_dec_new = nic_dec_new + 2*boil::BW +ii;
      int nj_dec_new = njc_dec_new + 2*boil::BW +jj;
      int nk_dec_new = nkc_dec_new + 2*boil::BW +kk;

      /*---------+
      |  output  |
      +---------*/
      /* allocate decomposed array */
      real *** val_new;
      alloc3d( &val_new, ni_dec_new, nj_dec_new, nk_dec_new );

      int ni_block_new = dims[0];
      int nj_block_new = dims[1];
      int nk_block_new = dims[2];
      cout<<"# ni_block_new, nj_block_new, nk_block_new= "
               <<ni_block_new<<" "<<nj_block_new<<" "<<nk_block_new<<"\n";

      /*----------+
      |  compose  |
      +----------*/
      for (int iproc = 0; iproc < np_new; iproc++) {

        int i_block = real(iproc + 0)/real(nj_block_new*nk_block_new) + 1;
        int j_block = real(iproc + 0 - (i_block-1)*(nj_block_new*nk_block_new))
                    / real(nk_block_new) + 1;
        int k_block = real(iproc + 1 - (i_block-1)*(nj_block_new*nk_block_new)
                          -(j_block-1)*nk_block_new); 
        cout<<"iproc= "<<iproc<<" "<<i_block<<" "<<j_block<<" "<<k_block<<"\n";

        /* copy q to val_new */
        int ist = (i_block-1) * nic_dec_new;
        int jst = (j_block-1) * njc_dec_new;
        int kst = (k_block-1) * nkc_dec_new;
        //cout<<"ist,jst,kst = "<<ist<<" "<<jst<<" "<<kst<<"\n";
        int imax=0; 
        int jmax=0; 
        int kmax=0; 
        for(int i=0; i<ni_dec_new; i++)
          for(int j=0; j<nj_dec_new; j++)
            for(int k=0; k<nk_dec_new; k++) {
              val_new[i][j][k] = q[ist+i][jst+j][kst+k];
              if(ist+i>imax) imax=ist+i;
              if(jst+j>jmax) jmax=jst+j;
              if(kst+k>kmax) kmax=kst+k;
            }
        //cout<<"imax= "<<imax<<" "<<jmax<<" "<<kmax<<"\n";

        /* file name */
        name = name_file(nm, ".bck2", ts, iproc);

        /* open a file */
        ofstream out(name.c_str(), ios::binary);

        /* write val_new */
        cout<<"# output: "<<name<<"\n";
        out.write(reinterpret_cast<char *> (&ni_dec_new), sizeof(int));
        out.write(reinterpret_cast<char *> (&nj_dec_new), sizeof(int));
        out.write(reinterpret_cast<char *> (&nk_dec_new), sizeof(int));
        out.write(reinterpret_cast<char *> (val_new[0][0]),
                                ni_dec_new*nj_dec_new*nk_dec_new*sizeof(real));
      }

      /*-------------+
      |  deallocate  |
      +-------------*/
      dealloc3d( &val_new );
    }

    dealloc3d( &q );

    return;
}

/******************************************************************************/
void vector_treat( const char *nm, const int nic, const int njc, const int nkc
             , const int ts, const int np, const int np_new){
  /*----------------------------------------------+
  |  Vector                                       |
  +----------------------------------------------*/
    int BW = boil::BW;
    int ni[] = { nic+2*BW+1, nic+2*BW  , nic+2*BW  };
    int nj[] = { njc+2*BW  , njc+2*BW+1, njc+2*BW  };
    int nk[] = { nkc+2*BW  , nkc+2*BW  , nkc+2*BW+1};
    real *** u;
    alloc3d( &u, ni[0], nj[0], nk[0] );
    real *** v;
    alloc3d( &v, ni[1], nj[1], nk[1] );
    real *** w;
    alloc3d( &w, ni[2], nj[2], nk[2] );

    for(int i=0; i<ni[0]; i++)
      for(int j=0; j<nj[0]; j++)
        for(int k=0; k<nk[0]; k++)
          u[i][j][k]=1.0e+300;

    for(int i=0; i<ni[1]; i++)
      for(int j=0; j<nj[1]; j++)
        for(int k=0; k<nk[1]; k++)
          v[i][j][k]=1.0e+300;

    for(int i=0; i<ni[2]; i++)
      for(int j=0; j<nj[2]; j++)
        for(int k=0; k<nk[2]; k++)
          w[i][j][k]=1.0e+300;

    string name;
    /* get decomposed size */
    name = name_file(nm, ".bck", ts, 0);
    ifstream in(name.c_str(), ios::binary);

    /* stop if file is not present */
    if( in.rdstate() != 0 ) {
      cout << "failed to open " << name << endl;
      cout << "exiting!" << endl;
      exit(0);
    }

    int ni_dec[3], nj_dec[3], nk_dec[3];
    for(int m=0; m<3; m++) {
      in.read(reinterpret_cast<char *> (&ni_dec[m]), sizeof(int));
      in.read(reinterpret_cast<char *> (&nj_dec[m]), sizeof(int));
      in.read(reinterpret_cast<char *> (&nk_dec[m]), sizeof(int));
      cout<<"ni_dec= "<<ni_dec[m]<<" "<<nj_dec[m]<<" "<<nk_dec[m]<<"\n";
    }

    /* allocate decomposed array */
    real *** u_dec, *** v_dec, *** w_dec;
    alloc3d( &u_dec, ni_dec[0], nj_dec[0], nk_dec[0] );
    alloc3d( &v_dec, ni_dec[1], nj_dec[1], nk_dec[1] );
    alloc3d( &w_dec, ni_dec[2], nj_dec[2], nk_dec[2] );

    int nic_dec=ni_dec[0]-2*BW-1;
    int njc_dec=nj_dec[0]-2*BW;
    int nkc_dec=nk_dec[0]-2*BW;

    cout<<"# nic, njc, nkc= "
        <<nic<<" "<<njc<<" "<<nkc<<"\n";
    cout<<"# nic_dec, njc_dec, nkc_dec= "
        <<nic_dec<<" "<<njc_dec<<" "<<nkc_dec<<"\n";

    int ni_block = nic/nic_dec;
    int nj_block = njc/njc_dec;
    int nk_block = nkc/nkc_dec;

    cout<<"# ni_block, nj_block, nk_block= "
        <<ni_block<<" "<<nj_block<<" "<<nk_block<<"\n";

    /*----------+
    |  compose  |
    +----------*/
    for (int iproc = 0; iproc < np; iproc++) {

      /* block ID */
      int i_block = real(iproc + 0)/real(nj_block*nk_block) + 1;
      int j_block = real(iproc + 0 - (i_block-1)*(nj_block*nk_block))
                  / real(nk_block) + 1;
      int k_block = real(iproc + 1 - (i_block-1)*(nj_block*nk_block)
                        -(j_block-1)*nk_block);
      //cout<<"iproc= "<<iproc<<" "<<i_block<<" "<<j_block<<" "<<k_block<<"\n";

      /* file name */
      name = name_file(nm, ".bck", ts, iproc);

      /* open a file */
      ifstream in(name.c_str(), ios::binary);

      /* stop if file is not present */
      if( in.rdstate() != 0 ) {
        cout << "failed to open " << name << endl;
        cout << "exiting!" << endl;
        exit(0);
      }

      int ni_tmp[3], nj_tmp[3], nk_tmp[3];
      for(int m=0; m<3; m++) {
        in.read(reinterpret_cast<char *> (&ni_tmp[m]), sizeof(int));
        in.read(reinterpret_cast<char *> (&nj_tmp[m]), sizeof(int));
        in.read(reinterpret_cast<char *> (&nk_tmp[m]), sizeof(int));
        //cout<<"ni_tmp= "<<ni_tmp[m]<<" "<<nj_tmp[m]<<" "<<nk_tmp[m]<<"\n";
      }
      in.read(reinterpret_cast<char *> (u_dec[0][0]),
            ni_dec[0]*nj_dec[0]*nk_dec[0] * sizeof(real));
      in.read(reinterpret_cast<char *> (v_dec[0][0]),
            ni_dec[1]*nj_dec[1]*nk_dec[1] * sizeof(real));
      in.read(reinterpret_cast<char *> (w_dec[0][0]),
            ni_dec[2]*nj_dec[2]*nk_dec[2] * sizeof(real));

      /* copy val to q */
      int ist = (i_block-1) * nic_dec;
      int jst = (j_block-1) * njc_dec;
      int kst = (k_block-1) * nkc_dec;
      //cout<<"ist,jst,kst = "<<ist<<" "<<jst<<" "<<kst<<"\n";

      int istl=boil::BW; int jstl=boil::BW; int kstl=boil::BW;
      if (i_block==1) istl=0;
      if (j_block==1) jstl=0;
      if (k_block==1) kstl=0;
      for(int i=istl; i<ni_dec[0]; i++)
        for(int j=jstl; j<nj_dec[0]; j++)
          for(int k=kstl; k<nk_dec[0]; k++) {
            u[ist+i][jst+j][kst+k] = u_dec[i][j][k];
          }

      istl=boil::BW; jstl=boil::BW; kstl=boil::BW;
      if (i_block==1) istl=0;
      if (j_block==1) jstl=0;
      if (k_block==1) kstl=0;
      for(int i=istl; i<ni_dec[1]; i++)
        for(int j=jstl; j<nj_dec[1]; j++)
          for(int k=kstl; k<nk_dec[1]; k++) {
            v[ist+i][jst+j][kst+k] = v_dec[i][j][k];
          }

      istl=boil::BW; jstl=boil::BW; kstl=boil::BW;
      if (i_block==1) istl=0;
      if (j_block==1) jstl=0;
      if (k_block==1) kstl=0;
      for(int i=istl; i<ni_dec[2]; i++)
        for(int j=jstl; j<nj_dec[2]; j++)
          for(int k=kstl; k<nk_dec[2]; k++) {
            w[ist+i][jst+j][kst+k] = w_dec[i][j][k];
          }
    }
    /*-------------+
    |  deallocate  |
    +-------------*/
    dealloc3d( &u_dec );
    dealloc3d( &v_dec );
    dealloc3d( &w_dec );

    /*------------------------+
    |  output to single-file  |
    +------------------------*/
    if (np_new==1) {
      /* file name */
      name = name_file(nm, ".bck2", ts, 0);

      /* open a file */
      ofstream out(name.c_str(), ios::binary);
      for(int m=0; m<3; m++) {
        out.write(reinterpret_cast<char *> (&ni[m]), sizeof(int));
        out.write(reinterpret_cast<char *> (&nj[m]), sizeof(int));
        out.write(reinterpret_cast<char *> (&nk[m]), sizeof(int));
      }
      out.write(reinterpret_cast<char *> (u[0][0]),
              ni[0]*nj[0]*nk[0] * sizeof(real));
      out.write(reinterpret_cast<char *> (v[0][0]),
              ni[1]*nj[1]*nk[1] * sizeof(real));
      out.write(reinterpret_cast<char *> (w[0][0]),
              ni[2]*nj[2]*nk[2] * sizeof(real));

    } else {

      /*------------+
      |  decompose  |
      +------------*/
      int * dims = new int[3];
      decomp(dims, np_new, nic, njc, nkc);
      cout<<"# dims[]= "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"\n";
      int ni_block_new = dims[0];
      int nj_block_new = dims[1];
      int nk_block_new = dims[2];
      cout<<"# ni_block_new, nj_block_new, nk_block_new= "
               <<ni_block_new<<" "<<nj_block_new<<" "<<nk_block_new<<"\n";

      int nic_dec_new=nic/dims[0];
      int njc_dec_new=njc/dims[1];
      int nkc_dec_new=nkc/dims[2];
      cout<<"# nic_dec_new= "<<nic_dec_new<<" "<<njc_dec_new<<" "
          <<nkc_dec_new<<"\n";
  
      int ni_dec_new[3], nj_dec_new[3], nk_dec_new[3];
      ni_dec_new[0] = nic_dec_new + 2*BW+1;
      nj_dec_new[0] = njc_dec_new + 2*BW;
      nk_dec_new[0] = nkc_dec_new + 2*BW;
  
      ni_dec_new[1] = nic_dec_new + 2*BW;
      nj_dec_new[1] = njc_dec_new + 2*BW+1;
      nk_dec_new[1] = nkc_dec_new + 2*BW;

      ni_dec_new[2] = nic_dec_new + 2*BW;
      nj_dec_new[2] = njc_dec_new + 2*BW;
      nk_dec_new[2] = nkc_dec_new + 2*BW+1;

      /*---------+
      |  output  |
      +---------*/
      /* allocate decomposed array */
      real *** u_dec_new, *** v_dec_new, *** w_dec_new;
      alloc3d( &u_dec_new, ni_dec_new[0], nj_dec_new[0], nk_dec_new[0] );
      alloc3d( &v_dec_new, ni_dec_new[1], nj_dec_new[1], nk_dec_new[1] );
      alloc3d( &w_dec_new, ni_dec_new[2], nj_dec_new[2], nk_dec_new[2] );
  
      /*------------+
      |  decompose  |
      +------------*/
      for (int iproc = 0; iproc < np_new; iproc++) {
  
        int i_block = real(iproc + 0)/real(nj_block_new*nk_block_new) + 1;
        int j_block = real(iproc + 0 - (i_block-1)*(nj_block_new*nk_block_new))
                    / real(nk_block_new) + 1;
        int k_block = real(iproc + 1 - (i_block-1)*(nj_block_new*nk_block_new)
                          -(j_block-1)*nk_block_new);
        //cout<<"iproc= "<<iproc<<" "<<i_block<<" "<<j_block<<" "<<k_block<<"\n";
  
        /* copy q to val_new */
        int ist = (i_block-1) * nic_dec_new;
        int jst = (j_block-1) * njc_dec_new;
        int kst = (k_block-1) * nkc_dec_new;
        //cout<<"ist= "<<ist<<" "<<jst<<" "<<kst<<"\n";
  
        for(int i=0; i<ni_dec_new[0]; i++)
          for(int j=0; j<nj_dec_new[0]; j++)
            for(int k=0; k<nk_dec_new[0]; k++)
              u_dec_new[i][j][k] = u[ist+i][jst+j][kst+k];
  
        for(int i=0; i<ni_dec_new[1]; i++)
          for(int j=0; j<nj_dec_new[1]; j++)
            for(int k=0; k<nk_dec_new[1]; k++)
              v_dec_new[i][j][k] = v[ist+i][jst+j][kst+k];
  
        for(int i=0; i<ni_dec_new[2]; i++)
          for(int j=0; j<nj_dec_new[2]; j++)
            for(int k=0; k<nk_dec_new[2]; k++){
              w_dec_new[i][j][k] = w[ist+i][jst+j][kst+k];
            }
  
        /* file name */
        name = name_file(nm, ".bck2", ts, iproc);
  
        /* open a file */
        ofstream out(name.c_str(), ios::binary);
  
        /* write val_new */
        cout<<"# output: "<<name<<"\n";
        for(int m=0; m<3; m++) {
          out.write(reinterpret_cast<char *> (&ni_dec_new[m]), sizeof(int));
          out.write(reinterpret_cast<char *> (&nj_dec_new[m]), sizeof(int));
          out.write(reinterpret_cast<char *> (&nk_dec_new[m]), sizeof(int));
        }
        out.write(reinterpret_cast<char *> (u_dec_new[0][0]),
                ni_dec_new[0]*nj_dec_new[0]*nk_dec_new[0] * sizeof(real));
        out.write(reinterpret_cast<char *> (v_dec_new[0][0]),
                ni_dec_new[1]*nj_dec_new[1]*nk_dec_new[1] * sizeof(real));
        out.write(reinterpret_cast<char *> (w_dec_new[0][0]),
                ni_dec_new[2]*nj_dec_new[2]*nk_dec_new[2] * sizeof(real));
      }
  
        /*-------------+
        |  deallocate  |
        +-------------*/
        dealloc3d( &u_dec_new );
        dealloc3d( &v_dec_new );
        dealloc3d( &w_dec_new );
    }

    dealloc3d( &u );
    dealloc3d( &v );
    dealloc3d( &w );

    return;
}
/******************************************************************************/
void nucleation( const char *nm, const int nic, const int njc, const int nkc
             , const int ts, const int np, const int np_new){

  string f_source = name_file(nm, ".bck", ts, 0);

  for (int iproc = 0; iproc < np_new; iproc++) {
    string f_dest = name_file(nm, ".bck2", ts, iproc);
    cout<<"# output: "<<f_dest<<"\n";
    ifstream in (f_source.c_str(), ios::binary);
    ofstream out(f_dest.c_str(),   ios::binary);

    if(in.is_open() && out.is_open()) {
      while(!in.eof()) {
        out.put(in.get());
      }
    }

    //Close both files
    in.close();
    out.close();
  }

}
/******************************************************************************/
void factor(int n, int * factor, int * number) {

  (*number) = 0;

  int m=n-1;
  for(;;) {
    if(n % m == 0) {
      factor[(*number)++] = n/m;
      n = m;
    }
    m--;
    if(m == 0) break;
    if(m == 1) {
      factor[(*number)++] = n;
      break;
    }
  }
}

/******************************************************************************/
void decomp(int * dis, const int np_new
          , const int nic, const int njc, const int nkc ) {

  /* get the resolution in each direction */
  int res[] = { nic , njc , nkc };

  /*---------------------------+
  |  initialize distributions  |
  +---------------------------*/
  for(int m=0; m<3; m++)
    dis[m] = 1;

  /*---------------------------------+
  |  factorize number of processors  |
  +---------------------------------*/
  int factors[64], nfact;
  factor( np_new, factors, &nfact);

  /*-------------------------------------------------+
  |  assign factor to coordinate directions (dis[])  |
  +-------------------------------------------------*/
  for(int f=0; f<nfact; f++) {

    int maxr = boil::maxi( res[0], res[1], res[2] );

    bool done = false;

    /* try the biggest one */
    for(int d=0; d<3; d++) {
      if( maxr == res[d] ) {
        if( res[d]%factors[f] == 0 ) {
          dis[d] *= factors[f];
          res[d] /= factors[f];
          done = true;
          break;
        }
      }
    }

    /* if it fails, go for the one which is dividable */
    if( !done ) {
      for(int d=0; d<3; d++) {
        if( res[d]%factors[f] == 0 ) {
          dis[d] *= factors[f];
          res[d] /= factors[f];
          done = true;
          break;
        }
      }
    }

    assert(done);
  }
}
