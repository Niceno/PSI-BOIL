const real LX =   8.0;
const real LY =   2.0;
const real LZ =   4.0;

/*----------+
|  grid(s)  |
+----------*/
Grid1D gx(Range<real>( 0.0, LX ), 
          64, Periodic::yes());

Grid1D gy(Range<real>( -LY/2.0,  LY/2.0), 
          Range<real>( LY/256.0, LY/256.0 ),
          64, Periodic::no());

Grid1D gz(Range<real>( 0.0, LZ ), 
          64, Periodic::yes());

/*---------+
|  domain  |
+---------*/
Domain d(gx, gy, gz);
