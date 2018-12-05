#include "finescalar.h"

void FineScalar::output(int ii, int jj, int kk) {
  boil::oout <<"FineScalar::otp@ "<< ii<<" "<<jj<<" "<<kk<<" "<<adens[ii][jj][kk]
     <<" "<<(*phi)[ii][jj][kk]<<" | "
     <<" "<<value(ii,jj,kk,w())<<" "<<value(ii,jj,kk,e())
     <<" "<<value(ii,jj,kk,s())<<" "<<value(ii,jj,kk,n())
     <<" "<<value(ii,jj,kk,b())<<" "<<value(ii,jj,kk,t())<<" |"
     <<" "<<value(ii,jj,kk,ws())<<" "<<value(ii,jj,kk,wn())
     <<" "<<value(ii,jj,kk,wb())<<" "<<value(ii,jj,kk,wt())<<" |"
     <<" "<<value(ii,jj,kk,es())<<" "<<value(ii,jj,kk,en())
     <<" "<<value(ii,jj,kk,eb())<<" "<<value(ii,jj,kk,et())<<" |"
     <<" "<<value(ii,jj,kk,sb())<<" "<<value(ii,jj,kk,st())
     <<" "<<value(ii,jj,kk,nb())<<" "<<value(ii,jj,kk,nt())<<" ||"
     <<" "<<value(ii,jj,kk,wsb())<<" "<<value(ii,jj,kk,wst())
     <<" "<<value(ii,jj,kk,wnb())<<" "<<value(ii,jj,kk,wnt())
     <<" "<<value(ii,jj,kk,esb())<<" "<<value(ii,jj,kk,est())
     <<" "<<value(ii,jj,kk,enb())<<" "<<value(ii,jj,kk,ent())
     <<boil::endl;
 
  return;
}
