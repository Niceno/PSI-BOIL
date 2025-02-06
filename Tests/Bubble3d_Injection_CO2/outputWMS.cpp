    /* output data at WMS */
    if( time.current_time() / period_WMS >= real(itWMS) ) {
      for(int iWMS = 0; iWMS<nWMS; iWMS++) {
        real cWMS[NX][NX];
        real tprWMS[NX][NX];
        real wWMS[NX][NX];
        for (int I = 0; I < NX; I++) {
          for (int J = 0; J < NX; J++) {
            cWMS[I][J]   = 0.0;
            tprWMS[I][J] = 0.0;
            wWMS[I][J]   = 0.0;
          }
        }
        real zcWMS = -10.0;  // cell cener coordinate
        for_vk(c,k) {
          if (c.zn(k)<=zWMS[iWMS] && zWMS[iWMS]<c.zn(k+1)) {
            zcWMS = c.zc(k);
            for_vij(c,i,j) {
              int I = d.global_I(i)-boil::BW;
              int J = d.global_J(j)-boil::BW;
              cWMS[I][J]   = c[i][j][k];
              tprWMS[I][J] = tpr[i][j][k];
              wWMS[I][J]   = uvw[Comp::w()][i][j][k];
            }
          }
        }
        boil::cart.max_real(&zcWMS);

        real c2D[NX*NX],tpr2D[NX*NX],w2D[NX*NX];
        for (int I = 0; I < NX*NX; I++) {
          c2D[I]   = 0.0;
          tpr2D[I] = 0.0;
          w2D[I]   = 0.0;
        }

        int Itmp = 0;
        for (int I = 0; I < NX; I++) {
          for (int J = 0; J < NX; J++) {
            c2D[Itmp]   = cWMS[I][J];
            tpr2D[Itmp] = tprWMS[I][J];
            w2D[Itmp]   = wWMS[I][J];
            Itmp++;
          }
        }
        boil::cart.sum_real_n(c2D,NX*NX);
        boil::cart.sum_real_n(tpr2D,NX*NX);
        boil::cart.sum_real_n(w2D,NX*NX);

        Itmp = 0;
        for (int I = 0; I < NX; I++) {
          for (int J = 0; J < NX; J++) {
            cWMS[I][J]   = c2D[Itmp];
            tprWMS[I][J] = tpr2D[Itmp];
            wWMS[I][J]   = w2D[Itmp];
            Itmp++;
          }
        }

        // output to the file data2D-itWMS.txt
        if (boil::cart.iam()==0) {
          std::string name0 = "WMS"+std::to_string(iWMS);
          std::string name = name_file(name0.c_str(), ".dat", itWMS);
          std::fstream output;
          output << std::setprecision(16);
          output.open(name, std::ios::out);
          output << "###time= " << time.current_time() << " frame= " << itWMS <<"\n";
          output << "VARIABLES= X Y Z vfl eps W\n";
          output << "ZONE I= "<<NX<<" ,J= "<<NX<<" ,K=1 DATAPACKING=POINT\n";

          for (int I = 0; I < NX; I++) {
            for (int J = 0; J < NX; J++) {
              output<<d.xc_global(I+boil::BW)<<" "<<d.yc_global(J+boil::BW)<<" "<<zcWMS
              <<" "<<cWMS[I][J]<<" "<<tprWMS[I][J]<<" "<<wWMS[I][J]<<"\n";
            }
          }
          output.close();
        }
      }
      // counter ++ 
      itWMS = int(time.current_time()/period_WMS) + 1;
    }

