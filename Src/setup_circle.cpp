void setup_circle(Scalar & c, const real radius, const real xcent, const real ycent) {  
#if 1
  // set theoretical value of single droplet for color function
  for_vijk(c,i,j,k){
    int inm,inp,jnm,jnp;
    if (c.xc(i)>xcent && c.yc(j)>ycent) {
      inm=i; inp=i+1; jnm=j; jnp=j+1;
    } else if (c.xc(i)>xcent && c.yc(j)<ycent) {
      inm=i; inp=i+1; jnm=j+1; jnp=j; 
    } else if (c.xc(i)<xcent && c.yc(j)>ycent) {
      inm=i+1; inp=i; jnm=j; jnp=j+1;
    } else {
      inm=i+1; inp=i; jnm=j+1; jnp=j;
    }

    // check vertex
    int imm=0, imp=0, ipm=0, ipp=0;

    real dtmp = sqrt(pow(c.xn(inm)-xcent,2.0) + pow(c.yn(jnm)-ycent,2.0));
    if (dtmp<=radius) imm=1;

    dtmp = sqrt(pow(c.xn(inm)-xcent,2.0) + pow(c.yn(jnp)-ycent,2.0));
    if (dtmp<=radius) imp=1;

    dtmp = sqrt(pow(c.xn(inp)-xcent,2.0) + pow(c.yn(jnm)-ycent,2.0));
    if (dtmp<=radius) ipm=1;

    dtmp = sqrt(pow(c.xn(inp)-xcent,2.0) + pow(c.yn(jnp)-ycent,2.0));
    if (dtmp<=radius) ipp=1;

    if (imm+imp+ipm+ipp==4) {

      // full of liquid
      c[i][j][k] = 1.0;

    } else if(imm+imp+ipm+ipp==0){

      // full of gas
      c[i][j][k] = 0.0;

    } else if(imm+imp+ipm+ipp==1) {

      // one vertex is in liquid
      real theta1 = asin(fabs(c.yn(jnm)-ycent)/radius);
      real theta2 = acos(fabs(c.xn(inm)-xcent)/radius);

      if (approx(theta2-theta1,0.0)) {

        //exception
        c[i][j][k]=0.0;
      } else {

        real area = radius*radius*(theta2-theta1)*0.5;
        real a = fabs(c.yn(jnm)-ycent)/(radius*cos(theta1));
        real yy = a*fabs(c.xn(inm)-xcent);
        area = area - (radius*sin(theta2)-yy)*fabs(c.xn(inm)-xcent)*0.5
                    - (fabs(c.yn(jnm)-ycent)-yy)
                     *(radius*cos(theta1)-fabs(c.xn(inm)-xcent))*0.5;
        c[i][j][k] = area/(fabs(c.dxc(inm))*fabs(c.dyc(jnm)));
        //std::cout<<area<<"\n";
        
      }
#if 1
    } else if (imm+imp+ipm+ipp==3) {

      // one vertex is in gas
      real theta1 = acos(fabs(c.xn(inp)-xcent)/radius);
      real theta2 = asin(fabs(c.yn(jnp)-ycent)/radius);

      real area = (fabs(c.xn(inp)-xcent)-radius*cos(theta2))*fabs(c.yn(jnp)-ycent)*0.5
                + (fabs(c.yn(jnp)-ycent)-radius*sin(theta1))*fabs(c.xn(inp)-xcent)*0.5;
      area -= radius*radius*(theta2-theta1)*0.5;
      area = c.dxc(i)*c.dyc(j)-area;
      c[i][j][k] = area/(c.dxc(i)*c.dyc(j));
      //std::cout<<c[i][j][k]<<"\n";

    } else if (imm+imp+ipm+ipp==2) {

      // two vertices are in gas
      if (ipm==1) {
        // imm and ipm are in liquid
        real theta1 = acos(fabs(c.xn(inp)-xcent)/radius);
        real theta2 = acos(fabs(c.xn(inm)-xcent)/radius);

        real area = radius*radius*(theta2-theta1)*0.5;
        real a1 = radius*sin(theta1)/fabs(c.xn(inp)-xcent);
        area += (fabs(c.xn(inp)-xcent)-fabs(c.yn(jnm)-ycent)/a1)
               *(radius*sin(theta1)-fabs(c.yn(jnm)-ycent))*0.5;
        area -= (fabs(c.yn(jnm)-ycent)/a1-fabs(c.xn(inm)-xcent))*fabs(c.yn(jnm)-ycent)*0.5;
        area -= (radius*sin(theta2)-fabs(c.yn(jnm)-ycent))*fabs(c.xn(inm)-xcent)*0.5;
        c[i][j][k] = area/(c.dxc(i)*c.dyc(j));
        //std::cout<<c[i][j][k]<<"\n";
      } else if (imp==1) {
        // imm and imp are in liquid
        real theta1 = asin(fabs(c.yn(jnm)-ycent)/radius);
        real theta2 = asin(fabs(c.yn(jnp)-ycent)/radius);

        real area = radius*radius*(theta2-theta1)*0.5;
        real a2 = fabs(c.yn(jnp)-ycent)/(radius*cos(theta2));
        area += (radius*cos(theta2)-fabs(c.xn(inm)-xcent))
               *(fabs(c.yn(jnp)-ycent)-a2*fabs(c.xn(inm)-xcent))*0.5;
        area -= (a2*fabs(c.xn(inm)-xcent)-fabs(c.yn(jnm)-ycent))*fabs(c.xn(inm)-xcent)*0.5;
        area -= (radius*cos(theta1)-fabs(c.xn(inm)-xcent))*fabs(c.yn(jnm)-ycent)*0.5;
        c[i][j][k] = area/(c.dxc(i)*c.dyc(j));
      } else {
        std::cout<<"Error!!!\n";
        exit(0);
      }
#endif
    }
  }
#else
  for_vijk(c,i,j,k) {
    real wsb_x = c.xc(i) - c.dxc(i)*0.5;
    real wst_x = c.xc(i) - c.dxc(i)*0.5;
    real wnb_x = c.xc(i) - c.dxc(i)*0.5;
    real wnt_x = c.xc(i) - c.dxc(i)*0.5;
    real esb_x = c.xc(i) + c.dxc(i)*0.5;
    real est_x = c.xc(i) + c.dxc(i)*0.5;
    real enb_x = c.xc(i) + c.dxc(i)*0.5;
    real ent_x = c.xc(i) + c.dxc(i)*0.5;

    real wsb_y = c.yc(j) - c.dyc(j)*0.5;
    real wst_y = c.yc(j) - c.dyc(j)*0.5;
    real wnb_y = c.yc(j) + c.dyc(j)*0.5;
    real wnt_y = c.yc(j) + c.dyc(j)*0.5;
    real esb_y = c.yc(j) - c.dyc(j)*0.5;
    real est_y = c.yc(j) - c.dyc(j)*0.5;
    real enb_y = c.yc(j) + c.dyc(j)*0.5;
    real ent_y = c.yc(j) + c.dyc(j)*0.5;

    real wsb_z = c.zc(k) - c.dzc(k)*0.5;
    real wst_z = c.zc(k) + c.dzc(k)*0.5;
    real wnb_z = c.zc(k) - c.dzc(k)*0.5;
    real wnt_z = c.zc(k) + c.dzc(k)*0.5;
    real esb_z = c.zc(k) - c.dzc(k)*0.5;
    real est_z = c.zc(k) + c.dzc(k)*0.5;
    real enb_z = c.zc(k) - c.dzc(k)*0.5;
    real ent_z = c.zc(k) + c.dzc(k)*0.5;

    real wsb_dist = sqrt(pow(wsb_x-xcent,2.0)+pow(wsb_y-ycent,2.0));
    real wst_dist = sqrt(pow(wst_x-xcent,2.0)+pow(wst_y-ycent,2.0));
    real wnb_dist = sqrt(pow(wnb_x-xcent,2.0)+pow(wnb_y-ycent,2.0));
    real wnt_dist = sqrt(pow(wnt_x-xcent,2.0)+pow(wnt_y-ycent,2.0));
    real esb_dist = sqrt(pow(esb_x-xcent,2.0)+pow(esb_y-ycent,2.0));
    real est_dist = sqrt(pow(est_x-xcent,2.0)+pow(est_y-ycent,2.0));
    real enb_dist = sqrt(pow(enb_x-xcent,2.0)+pow(enb_y-ycent,2.0));
    real ent_dist = sqrt(pow(ent_x-xcent,2.0)+pow(ent_y-ycent,2.0));

    if(wsb_dist<radius&&wst_dist<radius&&wnb_dist<radius&&wnt_dist<radius&&
       esb_dist<radius&&est_dist<radius&&enb_dist<radius&&ent_dist<radius) {
       c[i][j][k] = 1.0;
    } else if(wsb_dist<=radius||wst_dist<=radius||wnb_dist<=radius||wnt_dist<=radius||
              esb_dist<=radius||est_dist<=radius||enb_dist<=radius||ent_dist<=radius) {
       int mm=30;
       real x0=c.xn(i);
       real y0=c.yn(j);
       real z0=c.zn(k);
       real ddx=c.dxc(i)/real(mm);
       real ddy=c.dyc(j)/real(mm);
       real ddz=c.dzc(k)/real(mm);
       int itmp=0;
       for (int ii=0; ii<mm; ii++){
         for (int jj=0; jj<mm; jj++){
           for (int kk=0; kk<mm; kk++){
             real xxc=x0+0.5*ddx+real(ii)*ddx;
             real yyc=y0+0.5*ddy+real(jj)*ddy;
             real zzc=z0+0.5*ddz+real(kk)*ddz;
             real dist=sqrt(pow(xxc-xcent,2.0)
                           +pow(yyc-ycent,2.0));
             if (dist<radius){
               itmp=itmp+1;
             }
           }
         }
       }
       c[i][j][k]=real(itmp)/real(mm*mm*mm);
    }
  }
#endif

  return;
}
