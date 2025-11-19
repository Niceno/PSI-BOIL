    /* set nucleation sites */
    srand(0);
    real x_min_nsd = -half_heater+1e-3;
    real x_max_nsd =  half_heater-1e-3;
    real y_min_nsd = -half_heater+1e-3;
    real y_max_nsd =  half_heater-1e-3;
    real zns = rseed * cos(45.0/180.0*boil::pi);  // Please edit
    real t_act = 40.0;		                  // Please edit
    fstream opt_nsd;
    opt_nsd.open("site.txt", ios::out);
    opt_nsd << "X Y Tact Z ID NSD\n";

    for(int nns=1; nns<=30; nns++) {

      /* set location of site */
      /* create 8 candidates */
      real xx[8], yy[8], dnear[8];  // candidates
      for(int i=0; i<8; i++){
        xx[i] = x_min_nsd + (x_max_nsd - x_min_nsd)*(double)rand()/(double)RAND_MAX;
        yy[i] = y_min_nsd + (y_max_nsd - y_min_nsd)*(double)rand()/(double)RAND_MAX;
      }

      /* find farest xns and yns */
      real xns, yns;  // selected seed point
      if(nns==1) {
        xns = xx[0];
        yns = yy[0];
      } else {
        for(int i=0; i<8; i++){
          dnear[i]=boil::exa;
	  /* distance between i and existing sites */
          for(int nsd=0; nsd<nucl.size(); nsd++) {
            real xtmp = nucl.sites[nsd].x();
            real ytmp = nucl.sites[nsd].y();
            real dd=sqrt((xtmp-xx[i])*(xtmp-xx[i])+(ytmp-yy[i])*(ytmp-yy[i]));
            if(dd<dnear[i]){
              dnear[i]=dd;
            }
          }
	  /* distance between i and existing dummy sites */
          for(int nsd=0; nsd<nucl.dsize(); nsd++) {
            real xtmp = nucl.dsites[nsd].x();
            real ytmp = nucl.dsites[nsd].y();
            real dd=sqrt((xtmp-xx[i])*(xtmp-xx[i])+(ytmp-yy[i])*(ytmp-yy[i]));
            if(dd<dnear[i]){
              dnear[i]=dd;
            }
          }
        }
        real dist=dnear[0];
        int ii=0;
        for(int i=1; i<8; i++){
          if(dist<dnear[i]) {
            dist=dnear[i];
            ii=i;
          }
        }
        //std::cout<<"dist,ii= "<<dist<<" "<<ii<<"\n";
        xns = xx[ii];
        yns = yy[ii];
      }

      nucl.add(Site( xns,  yns, zns, t_act, -0.002));
      boil::oout<<"main:nucl.add= "<<nns-1<<" "<<xns<<" "<<yns<<" "<<t_act<<"\n";
      opt_nsd << xns <<" "<<yns<<" "<<t_act<<" "<<zns<<" "<<nns<<"\n";
    }

    opt_nsd.close();

