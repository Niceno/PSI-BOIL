const real LX =   8.0;
const real LY =   2.0;
const real LZ =   4.0;

/*----------+
|  grid(s)  |
+----------*/
Grid1D gx(Range<real>( 0.0, LX ), 
          128, Periodic::yes());

Grid1D gy(Range<real>( -LY/2.0,  LY/2.0), 
          Range<real>( LY/512.0, LY/512.0 ),
          128, Periodic::no());

Grid1D gz(Range<real>( 0.0, LZ ), 
          128, Periodic::yes());

/*---------+
|  domain  |
+---------*/
Domain d(gx, gy, gz);
