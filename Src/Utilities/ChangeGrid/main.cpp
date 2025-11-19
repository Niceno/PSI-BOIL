#include "changegrid.h"
using namespace std;

void load_grid_nijk(ifstream & in, int & ni, int & nj, int & nk);
void load_grid_xyz(ifstream & in, int & ni, int & nj, int & nk,
                  double * xn, double * yn, double * zn);
void factor(int n, int * factor, int * number);
void decomp(int * id, const int i1, const int i2, const int i3, const int i4 );
void scalar( const char *nm,
             const int * ic2ic, const int * jc2jc, const int *kc2kc,
             const int ts, const int np, const int np_new,
             const int ii, const int jj, const int kk);
void vector_treat( const char *nm,
               const int * ic2ic, const int * jc2jc, const int *kc2kc,
               const int * ic2ic_stg, const int * jc2jci_stg, const int *kc2kc_stg,
               const int ts, const int np, const int np_new);
void nucleation( const char *nm, const int ts, const int np_new);

int np, np_new; // number of process original and new
int ni_org,nj_org,nk_org; // ni, nj, nk of original grid
int nic_org, njc_org, nkc_org; // nic, njc, nkc of original grid
int ni_new,nj_new,nk_new; // ni, nj, nk of new grid
int nic_new,njc_new,nkc_new; // nic, njc, nkc of new grid
int ni_stg_org,nj_stg_org,nk_stg_org; // ni, nj, nk of original grid
int ni_stg_new,nj_stg_new,nk_stg_new; // ni, nj, nk of new grid
/******************************************************************************/
int main() {

   cout<<"### Change grid resolution ###\n";
   cout<<"# Input1: Original grid, new grid\n";
   cout<<"# Input2: bck files on original grid\n";
   cout<<"# Output: bck files on new grid\n";

   cout<<"### Reading input deck. ###\n";

  ifstream input_deck("input.txt");
  istream * inp = &cin;

  if(input_deck.is_open()) {
    inp  = &input_deck;
  } else {
    cout<<"# Input file *input.txt* not found. "
             <<"Defaulting to standard input.\n";
  }

   /****************************************************************************
   * Read original and new grid                                                *
   ****************************************************************************/
   /******************
   *  Original grid  *
   ******************/
   /* original file name */
   cout<<"Enter original grid file name\n";
   string fn_org;
   *inp >> fn_org;
   cout<<"Original grid file name= "<<fn_org<<"\n";

   /* open original grid */
   ifstream ifs_org(fn_org,ios::binary);

   if(ifs_org.rdstate() != 0) {
     cout << "failed to open " << fn_org << std::endl;
     cout << "exiting!" << std::endl;
     exit(0);
   }
   /* read ni,nj,nk of original grid */
   load_grid_nijk(ifs_org,ni_org,nj_org,nk_org);
   std::cout<<"orgiginal ni,nj,nk= "<<ni_org<<" "<<nj_org<<" "<<nk_org<<"\n";
   nic_org = ni_org - 2*boil::BW;
   njc_org = nj_org - 2*boil::BW;
   nkc_org = nk_org - 2*boil::BW;
   std::cout<<"orgiginal nic,njc,nkc= "<<nic_org<<" "<<njc_org<<" "<<nkc_org<<"\n";


   /* allocate original grid */
   double xn_org[ni_org];
   double yn_org[nj_org];
   double zn_org[nk_org];

   /* read xn,yn,zn of original grid */
   ifs_org.seekg(0, ios_base::beg);
   load_grid_xyz(ifs_org,ni_org,nj_org,nk_org,xn_org,yn_org,zn_org);

   /* allocate new grid */
   ni_stg_org = ni_org+1;
   nj_stg_org = nj_org+1;
   nk_stg_org = nk_org+1;
   double xn_stg_org[ni_stg_org];
   double yn_stg_org[nj_stg_org];
   double zn_stg_org[nk_stg_org];

   /* read xn_stg,yn_stg,zn_stg of new grid */
   ifs_org.seekg(0, ios_base::beg);
   load_grid_xyz(ifs_org,ni_stg_org,nj_stg_org,nk_stg_org,
                 xn_stg_org,yn_stg_org,zn_stg_org);

   /* close original grid */
   ifs_org.close();


   /*************
   *  New grid  *
   *************/
   /* new file name */
   cout<<"Enter new grid file name\n";
   string fn_new;
   *inp >> fn_new;
   cout<<"New grid file name= "<<fn_new<<"\n";

   /* open new grid */
   ifstream ifs_new(fn_new,ios::binary);

   if(ifs_new.rdstate() != 0) {
     cout << "failed to open " << fn_new << std::endl;
     cout << "exiting!" << std::endl;
     exit(0);
   }
   /* read ni,nj,nk of new grid */
   load_grid_nijk(ifs_new,ni_new,nj_new,nk_new);
   std::cout<<"new ni,nj,nk= "<<ni_new<<" "<<nj_new<<" "<<nk_new<<"\n";

   nic_new = ni_new - 2*boil::BW;
   njc_new = nj_new - 2*boil::BW;
   nkc_new = nk_new - 2*boil::BW;
   std::cout<<"new nic,njc,nkc= "<<nic_new<<" "<<njc_new<<" "<<nkc_new<<"\n";
   
   /* allocate new grid */
   double xn_new[ni_new];
   double yn_new[nj_new];
   double zn_new[nk_new];
   
   /* read xn,yn,zn of new grid */
   ifs_new.seekg(0, ios_base::beg);
   load_grid_xyz(ifs_new,ni_new,nj_new,nk_new,xn_new,yn_new,zn_new);

   /* allocate new grid */
   ni_stg_new = ni_new+1;
   nj_stg_new = nj_new+1;
   nk_stg_new = nk_new+1;
   double xn_stg_new[ni_stg_new];
   double yn_stg_new[nj_stg_new];
   double zn_stg_new[nk_stg_new];
   
   /* read xn_stg,yn_stg,zn_stg of new grid */
   ifs_new.seekg(0, ios_base::beg);
   load_grid_xyz(ifs_new,ni_stg_new,nj_stg_new,nk_stg_new,
                 xn_stg_new,yn_stg_new,zn_stg_new);

   /* close new grid */
   ifs_new.close();

   /****************************************************************************
   * Define pointer from new-cell-index to original-cell-index                 *
   ****************************************************************************/
   /*-------------------------+
   | ic2ic[inew] returns iorg |
   +-------------------------*/
   int ic2ic[ni_new];
   for(int inew=0; inew<ni_new-1; inew++) {
     double xc_new = 0.5*(xn_new[inew] + xn_new[inew+1]);
     if(abs(xn_new[inew]-xn_new[inew+1])<1e-16) {
       if (inew <3 ) {
         ic2ic[inew]=inew;
         continue;
       } else if (inew>3) {
         ic2ic[inew]=ni_org-(ni_new-inew);
	 continue;
       }
     }
     for(int iorg=0; iorg<ni_org-1; iorg++) {
       if(xn_org[iorg]<xc_new && xc_new<=xn_org[iorg+1]) {
         ic2ic[inew]=iorg;
	 continue;
       }
     }
   }
   ic2ic[ni_new-1]=ni_org-1;
   //for(int inew=0; inew<ni_new; inew++) {
   //  cout<<"iorg= "<<ic2ic[inew]<<" inew= "<<inew<<"\n";
   //}
   //exit(0);

   /*-------------------------+
   | jc2jc[jnew] returns jorg |
   +-------------------------*/
   int jc2jc[nj_new];
   for(int jnew=0; jnew<nj_new-1; jnew++) {
     double yc_new = 0.5*(yn_new[jnew] + yn_new[jnew+1]);
     if(abs(yn_new[jnew]-yn_new[jnew+1])<1e-16) {
       if (jnew <3 ) {
         jc2jc[jnew]=jnew;
         continue;
       } else if (jnew>3) {
         jc2jc[jnew]=nj_org-(nj_new-jnew);
         continue;
       }
     }
     for(int jorg=0; jorg<nj_org-1; jorg++) {
       if(yn_org[jorg]<yc_new && yc_new<=yn_org[jorg+1]) {
         jc2jc[jnew]=jorg;
         continue;
       }
     }
   }
   jc2jc[nj_new-1]=nj_org-1;
   //for(int jnew=0; jnew<nj_new-1; jnew++) {
   //  int jorg=jc2jc[jnew];
   //  cout<<"jorg= "<<jorg<<" jnew= "<<jnew<<" "<<yn_org[jorg]
   //  <<" "<<0.5*(yn_new[jnew]+yn_new[jnew+1])<<" "<<yn_org[jorg+1]<<"\n";
   //}

   /*-------------------------+
   | kc2kc[knew] returns korg |
   +-------------------------*/
   int kc2kc[nk_new];
   for(int knew=0; knew<nk_new-1; knew++) {
     double zc_new = 0.5*(zn_new[knew] + zn_new[knew+1]);
     if(abs(zn_new[knew]-zn_new[knew+1])<1e-16) {
       if (knew <3 ) {
         kc2kc[knew]=knew;
         continue;
       } else if (knew>3) {
         kc2kc[knew]=nk_org-(nk_new-knew);
         continue;
       }
     }
     for(int korg=0; korg<nk_org-1; korg++) {
       if(zn_org[korg]<zc_new && zc_new<=zn_org[korg+1]) {
         kc2kc[knew]=korg;
         continue;
       }
     }
   }
   kc2kc[nk_new-1]=nk_org-1;
   //for(int knew=0; knew<nk_new-1; knew++) {
   //  int korg=kc2kc[knew];
   //  cout<<"korg= "<<korg<<" knew= "<<knew<<" "<<zn_org[korg]
   //  <<" "<<0.5*(zn_new[knew]+zn_new[knew+1])<<" "<<zn_org[korg+1]<<"\n";
   //}

   /*-----------------------------+
   | ic2ic_stg[inew] returns iorg |
   +-----------------------------*/
   int ic2ic_stg[ni_stg_new];
   for(int inew_stg=0; inew_stg<ni_stg_new-1; inew_stg++) {
     double xc_stg_new = 0.5*(xn_stg_new[inew_stg] + xn_stg_new[inew_stg+1]);
     if(abs(xn_stg_new[inew_stg]-xn_stg_new[inew_stg+1])<1e-16) {
       if (inew_stg <3 ) {
         ic2ic_stg[inew_stg]=inew_stg;
         continue;
       } else if (inew_stg>3) {
         ic2ic_stg[inew_stg]=ni_stg_org-(ni_stg_new-inew_stg);
         continue;
       }
     }
     for(int iorg_stg=0; iorg_stg<ni_stg_org-1; iorg_stg++) {
       if(xn_stg_org[iorg_stg]<xc_stg_new && xc_stg_new<=xn_stg_org[iorg_stg+1]) {
         ic2ic_stg[inew_stg]=iorg_stg;
         continue;
       }
     }
   }
   ic2ic_stg[ni_stg_new-1]=ni_stg_org-1;
   //for(int inew=0; inew<ni_new; inew++) {
   //  cout<<"iorg= "<<ic2ic[inew]<<" inew= "<<inew<<"\n";
   //}
   //exit(0);

   /*---------------------------------+
   | jc2jc_stg[jnew_stg] returns jorg |
   +---------------------------------*/
   int jc2jc_stg[nj_stg_new];
   for(int jnew_stg=0; jnew_stg<nj_stg_new-1; jnew_stg++) {
     double yc_stg_new = 0.5*(yn_stg_new[jnew_stg] + yn_stg_new[jnew_stg+1]);
     if(abs(yn_stg_new[jnew_stg]-yn_stg_new[jnew_stg+1])<1e-16) {
       if (jnew_stg <3 ) {
         jc2jc_stg[jnew_stg]=jnew_stg;
         continue;
       } else if (jnew_stg>3) {
         jc2jc_stg[jnew_stg]=nj_stg_org-(nj_stg_new-jnew_stg);
         continue;
       }
     }
     for(int jorg_stg=0; jorg_stg<nj_stg_org-1; jorg_stg++) {
       if(yn_stg_org[jorg_stg]<yc_stg_new && yc_stg_new<=yn_stg_org[jorg_stg+1]) {
         jc2jc_stg[jnew_stg]=jorg_stg;
         continue;
       }
     }
   }
   jc2jc_stg[nj_stg_new-1]=nj_stg_org-1;
   //for(int jnew=0; jnew<nj_new-1; jnew++) {
   //  int jorg=jc2jc[jnew];
   //  cout<<"jorg= "<<jorg<<" jnew= "<<jnew<<" "<<yn_org[jorg]
   //  <<" "<<0.5*(yn_new[jnew]+yn_new[jnew+1])<<" "<<yn_org[jorg+1]<<"\n";
   //}

   /*-------------------------+
   | kc2kc_stg[knew_stg] returns korg |
   +-------------------------*/
   int kc2kc_stg[nk_stg_new];
   for(int knew_stg=0; knew_stg<nk_stg_new-1; knew_stg++) {
     double zc_stg_new = 0.5*(zn_stg_new[knew_stg] + zn_stg_new[knew_stg+1]);
     if(abs(zn_stg_new[knew_stg]-zn_stg_new[knew_stg+1])<1e-16) {
       if (knew_stg <3 ) {
         kc2kc_stg[knew_stg]=knew_stg;
         continue;
       } else if (knew_stg>3) {
         kc2kc_stg[knew_stg]=nk_stg_org-(nk_stg_new-knew_stg);
         continue;
       }
     }
     for(int korg_stg=0; korg_stg<nk_stg_org-1; korg_stg++) {
       if(zn_stg_org[korg_stg]<zc_stg_new && zc_stg_new<=zn_stg_org[korg_stg+1]) {
         kc2kc_stg[knew_stg]=korg_stg;
         continue;
       }
     }
   }
   kc2kc_stg[nk_stg_new-1]=nk_stg_org-1;
#if 0
   for(int knew_stg=0; knew_stg<nk_stg_new-1; knew_stg++) {
     int korg_stg=kc2kc_stg[knew_stg];
     cout<<"korg_stg= "<<korg_stg<<" knew_stg= "<<knew_stg<<" "<<zn_stg_org[korg_stg]
     <<" "<<0.5*(zn_stg_new[knew_stg]+zn_stg_new[knew_stg+1])<<" "
     <<zn_stg_org[korg_stg+1]<<"\n";
   }
#endif

   /****************************************************************************
   * Change grid                                                               *
   ****************************************************************************/
    cout<<"# Enter current np (No. processors).\n";
    *inp>>np;
    cout<<"# np= "<<np<<"\n";

    cout<<"# Enter new np (No. processors).\n";
    *inp>>np_new;
    cout<<"# np_new= "<<np_new<<"\n";

   int ts;
   cout<<"# Enter time step.\n";
   *inp >> ts;
   cout<<"# time step= "<<ts<<"\n";

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
     *inp >> itype;
     cout<<"# file type = "<<itype<<"\n";
     cout<<"### Input deck built successfully. ###\n";

     string str_tmp;

     switch(itype) {

       case 0 :
         scalar(nm, ic2ic, jc2jc, kc2kc, ts, np, np_new,0,0,0);
         break;
       case 1 :
         vector_treat(nm, ic2ic, jc2jc, kc2kc, ic2ic_stg, jc2jc_stg, kc2kc_stg,
                      ts, np, np_new);
         break;
       case 2 :
         /* phi: cell center */
         str_tmp = nm;
         str_tmp += "-phi";
         scalar( nm, ic2ic, jc2jc, kc2kc, ts, np, np_new,0,0,0 );
         //scalar ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new,0,0,0 );

         /* clr: cell center */
         str_tmp = nm;
         str_tmp += "-clr";
         //scalar( nm, ic2ic, jc2jc, kc2kc, ts, np, np_new,0,0,0 );
         //scalar ( str_tmp.c_str(), ni_org, njc, nkc, ts, np, np_new,0,0,0 );
         /* f: node (vertex) */
         str_tmp = nm;
         str_tmp += "-f";
         //scalar ( str_tmp.c_str(), ni_org, njc, nkc, ts, np, np_new,1,1,1 );
         /* sigx: edge-x */
         str_tmp = nm;
         str_tmp += "-sigx";
         //scalar ( str_tmp.c_str(), ni_org, njc, nkc, ts, np, np_new,0,1,1 );
         /* sigx: edge-y */
         str_tmp = nm;
         str_tmp += "-sigy";
         //scalar ( str_tmp.c_str(), ni_org, njc, nkc, ts, np, np_new,1,0,1 );
         /* sigx: edge-y */
         str_tmp = nm;
         str_tmp += "-sigz";
         //scalar ( str_tmp.c_str(), ni_org, njc, nkc, ts, np, np_new,1,1,0 );
         /* sxyz: vecto */
         str_tmp = nm;
         str_tmp += "-face";
         //vector_treat ( str_tmp.c_str(), nic, njc, nkc, ts, np, np_new );
         break;
       case 3 :
         nucleation ( nm, ts, np_new );
         break;
       default :
         cout<<" Unrecognized file type. Exiting."<<endl;
         exit(0);
     }

    int ians;
    cout<<"# continue? 1(yes) or 0(no).\n";
    //cin>>ians;
    *inp >> ians;
    if (ians!=1) break;
  }


}
/******************************************************************************/
void load_grid_nijk(ifstream & in, int & ni, int & nj, int & nk) {
  in.read(reinterpret_cast<char *> (&ni), sizeof(int));
  in.read(reinterpret_cast<char *> (&nj), sizeof(int));
  in.read(reinterpret_cast<char *> (&nk), sizeof(int));
}

/******************************************************************************/
void load_grid_xyz( ifstream & in, int & ni, int & nj, int & nk,
                   double * xn, double * yn, double * zn){
  in.read(reinterpret_cast<char *> (&ni), sizeof(int));
  in.read(reinterpret_cast<char *> (&nj), sizeof(int));
  in.read(reinterpret_cast<char *> (&nk), sizeof(int));

  double tmp;
  for(int i=0; i<ni; i++){
    in.read(reinterpret_cast<char *> (&tmp),sizeof(double));
    xn[i] = tmp;
    //cout<<"i,xn= "<<i<<" "<<xn[i]<<"\n";
  }
  for(int j=0; j<nj; j++){
    in.read(reinterpret_cast<char *> (&tmp),sizeof(double));
    yn[j] = tmp;
    //cout<<"j,yn= "<<j<<" "<<yn[i]<<"\n";
  }
  for(int k=0; k<nk; k++){
    in.read(reinterpret_cast<char *> (&tmp),sizeof(double));
    zn[k] = tmp;
    //cout<<"k,zn= "<<k<<" "<<zn[k]<<"\n";
  }
}

/*****************************************************************************/
void scalar( const char *nm, //const int ni, const int nj, const int nk, 
            // const int ni_new, const int nj_new, const int nk_new,
             const int * ic2ic, const int * jc2jc, const int *kc2kc,
             const int ts, const int np, const int np_new, 
             const int ii, const int jj, const int kk){
    /* allocate */
    //int ni = nic+2*boil::BW;
    //int nj = njc+2*boil::BW;
    //int nk = nkc+2*boil::BW;
    const int ni = ni_org;
    const int nj = nj_org;
    const int nk = nk_org;
    const int nic = nic_org;
    const int njc = njc_org;
    const int nkc = nkc_org;
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
    cout<<"# nic, njc, nkc= "<<nic<<" "<<njc<<" "<<nkc<<"\n";
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

    /*-------------------------------+
    |  create variable for new grid  |
    +-------------------------------*/
    cout<<"create variable for new grid\n";
    real *** q_new;
    alloc3d( &q_new, ni_new+ii, nj_new+jj, nk_new+kk );
    for(int i=0; i<ni_new; i++)
      for(int j=0; j<nj_new; j++)
        for(int k=0; k<nk_new; k++) 
          q_new[i][j][k]=1.0e+300;

    cout<<"copy q to q_new\n";
    for(int i=0; i<ni_new+ii; i++)
      for(int j=0; j<nj_new+jj; j++)
        for(int k=0; k<nk_new+kk; k++){
          int i_old = ic2ic[i];
          int j_old = jc2jc[j];
          int k_old = kc2kc[k];
	  //cout<<i<<" "<<j<<" "<<k<<" "<<i_old<<" "<<j_old<<" "<<k_old<<"\n";
	  q_new[i][j][k] = q[i_old][j_old][k_old];
	  if (abs(q_new[i][j][k])>1e+299) {
            cout<<"i,j,k,q_new= "<<i<<" "<<j<<" "<<k<<" "<<q_new[i][j][k]<<"\n";
	    exit(0);
	  }
	}

    /*----------------------------+
    |  output bck to single-file  |
    +----------------------------*/
    if (np_new==1) {
      /* file name */
      name = name_file(nm, ".bck2", ts, 0);

      int ni_out=ni_new+ii;
      int nj_out=nj_new+jj;
      int nk_out=nk_new+kk;

      cout<<"output:np,ni,nj,nk= "<<np_new<<" "<<ni_out<<" "<<nj_out<<" "<<nk_out<<"\n";

      /* open a file */
      ofstream out(name.c_str(), ios::binary);

      out.write(reinterpret_cast<const char *> (&ni_out), sizeof(int));
      out.write(reinterpret_cast<const char *> (&nj_out), sizeof(int));
      out.write(reinterpret_cast<const char *> (&nk_out), sizeof(int));
      out.write(reinterpret_cast<const char *> (q_new[0][0]),
              ni_out*nj_out*nk_out*sizeof(real));

    } else {

      /*------------+
      |  decompose  |
      +------------*/
      int * dims = new int[3];
      //int nic_new = ni_new - 2*boil::BW;
      //int njc_new = nj_new - 2*boil::BW;
      //int nkc_new = nk_new - 2*boil::BW;
      decomp(dims, np_new, nic_new, njc_new, nkc_new);
      cout<<"# dims[]= "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"\n";
      int nic_dec_new=nic_new/dims[0];
      int njc_dec_new=njc_new/dims[1];
      int nkc_dec_new=nkc_new/dims[2];
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

        /* copy q_new to val_new */
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
              val_new[i][j][k] = q_new[ist+i][jst+j][kst+k];
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
    dealloc3d( &q_new );

    return;
}

/******************************************************************************/
void vector_treat( const char *nm,
               const int * ic2ic, const int * jc2jc, const int *kc2kc,
               const int * ic2ic_stg, const int * jc2jc_stg, const int *kc2kc_stg,
               const int ts, const int np, const int np_new){
  /*----------------------------------------------+
  |  Vector                                       |
  +----------------------------------------------*/
    int BW = boil::BW;
    //int ni[] = { nic+2*BW+1, nic+2*BW  , nic+2*BW  };
    //int nj[] = { njc+2*BW  , njc+2*BW+1, njc+2*BW  };
    //int nk[] = { nkc+2*BW  , nkc+2*BW  , nkc+2*BW+1};
    int ni[] = { ni_org+1, ni_org, ni_org  };
    int nj[] = { nj_org, nj_org+1, nj_org  };
    int nk[] = { nk_org, nk_org  , nk_org+1};
    const int nic = nic_org;
    const int njc = njc_org;
    const int nkc = nkc_org;
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

    /*-------------------------------+
    |  create variable for new grid  |
    +-------------------------------*/
    cout<<"create variable for new grid\n";
    int ni_nw[] = { ni_new+1, ni_new, ni_new  };
    int nj_nw[] = { nj_new, nj_new+1, nj_new  };
    int nk_nw[] = { nk_new, nk_new  , nk_new+1};
    real *** u_new;
    alloc3d( &u_new, ni_nw[0], nj_nw[0], nk_nw[0] );
    real *** v_new;
    alloc3d( &v_new, ni_nw[1], nj_nw[1], nk_nw[1] );
    real *** w_new;
    alloc3d( &w_new, ni_nw[2], nj_nw[2], nk_nw[2] );

    for(int i=0; i<ni_nw[0]; i++)
      for(int j=0; j<nj_nw[0]; j++)
        for(int k=0; k<nk_nw[0]; k++)
          u_new[i][j][k]=1.0e+300;

    for(int i=0; i<ni_nw[1]; i++)
      for(int j=0; j<nj_nw[1]; j++)
        for(int k=0; k<nk_nw[1]; k++)
          v_new[i][j][k]=1.0e+300;

    for(int i=0; i<ni_nw[2]; i++)
      for(int j=0; j<nj_nw[2]; j++)
        for(int k=0; k<nk_nw[2]; k++)
          w_new[i][j][k]=1.0e+300;

    cout<<"copy u to u_new\n";
    for(int i=0; i<ni_nw[0]; i++)
      for(int j=0; j<nj_nw[0]; j++)
        for(int k=0; k<nk_nw[0]; k++){
          int i_old = ic2ic_stg[i];
          int j_old = jc2jc[j];
          int k_old = kc2kc[k];
          //cout<<i<<" "<<j<<" "<<k<<" "<<i_old<<" "<<j_old<<" "<<k_old<<"\n";
          u_new[i][j][k] = u[i_old][j_old][k_old];
        }

    cout<<"copy v to v_new\n";
    for(int i=0; i<ni_nw[1]; i++)
      for(int j=0; j<nj_nw[1]; j++)
        for(int k=0; k<nk_nw[1]; k++){
          int i_old = ic2ic[i];
          int j_old = jc2jc_stg[j];
          int k_old = kc2kc[k];
          //cout<<i<<" "<<j<<" "<<k<<" "<<i_old<<" "<<j_old<<" "<<k_old<<"\n";
          v_new[i][j][k] = v[i_old][j_old][k_old];
        }

    cout<<"copy w to w_new\n";
    for(int i=0; i<ni_nw[2]; i++)
      for(int j=0; j<nj_nw[2]; j++)
        for(int k=0; k<nk_nw[2]; k++){
          int i_old = ic2ic[i];
          int j_old = jc2jc[j];
          int k_old = kc2kc_stg[k];
          //cout<<i<<" "<<j<<" "<<k<<" "<<i_old<<" "<<j_old<<" "<<k_old<<"\n";
          w_new[i][j][k] = w[i_old][j_old][k_old];
        }

    /*------------------------+
    |  output to single-file  |
    +------------------------*/
    if (np_new==1) {
      /* file name */
      name = name_file(nm, ".bck2", ts, 0);

      /* open a file */
      ofstream out(name.c_str(), ios::binary);
      for(int m=0; m<3; m++) {
        out.write(reinterpret_cast<char *> (&ni_nw[m]), sizeof(int));
        out.write(reinterpret_cast<char *> (&nj_nw[m]), sizeof(int));
        out.write(reinterpret_cast<char *> (&nk_nw[m]), sizeof(int));
      }
      out.write(reinterpret_cast<char *> (u_new[0][0]),
              ni_nw[0]*nj_nw[0]*nk_nw[0] * sizeof(real));
      out.write(reinterpret_cast<char *> (v_new[0][0]),
              ni_nw[1]*nj_nw[1]*nk_nw[1] * sizeof(real));
      out.write(reinterpret_cast<char *> (w_new[0][0]),
              ni_nw[2]*nj_nw[2]*nk_nw[2] * sizeof(real));

    } else {

      /*------------+
      |  decompose  |
      +------------*/
      int * dims = new int[3];
      decomp(dims, np_new, nic_new, njc_new, nkc_new);
      cout<<"# dims[]= "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"\n";
      int ni_block_new = dims[0];
      int nj_block_new = dims[1];
      int nk_block_new = dims[2];
      cout<<"# ni_block_new, nj_block_new, nk_block_new= "
               <<ni_block_new<<" "<<nj_block_new<<" "<<nk_block_new<<"\n";

      int nic_dec_new=nic_new/dims[0];
      int njc_dec_new=njc_new/dims[1];
      int nkc_dec_new=nkc_new/dims[2];
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
              u_dec_new[i][j][k] = u_new[ist+i][jst+j][kst+k];
  
        for(int i=0; i<ni_dec_new[1]; i++)
          for(int j=0; j<nj_dec_new[1]; j++)
            for(int k=0; k<nk_dec_new[1]; k++)
              v_dec_new[i][j][k] = v_new[ist+i][jst+j][kst+k];
  
        for(int i=0; i<ni_dec_new[2]; i++)
          for(int j=0; j<nj_dec_new[2]; j++)
            for(int k=0; k<nk_dec_new[2]; k++){
              w_dec_new[i][j][k] = w_new[ist+i][jst+j][kst+k];
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
    dealloc3d( &u_new );
    dealloc3d( &v_new );
    dealloc3d( &w_new );


    return;
}
/******************************************************************************/
void nucleation( const char *nm, const int ts, const int np_new){

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
