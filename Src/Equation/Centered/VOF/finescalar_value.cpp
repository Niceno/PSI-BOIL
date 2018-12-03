#include "finescalar.h"

/* interface to components */
real & FineScalar::value(int i, int j, int k, int dir) {
  if(dir == wsb()) return nodeval[i  ][j  ][k  ];
  if(dir == wst()) return nodeval[i  ][j  ][k+1];
  if(dir == wnb()) return nodeval[i  ][j+1][k  ];
  if(dir == wnt()) return nodeval[i  ][j+1][k+1];
  if(dir == esb()) return nodeval[i+1][j  ][k  ];
  if(dir == est()) return nodeval[i+1][j  ][k+1];
  if(dir == enb()) return nodeval[i+1][j+1][k  ];
  if(dir == ent()) return nodeval[i+1][j+1][k+1];

  if(dir == ws()) return edgez[i  ][j  ][k  ];
  if(dir == wn()) return edgez[i  ][j+1][k  ];
  if(dir == es()) return edgez[i+1][j  ][k  ];
  if(dir == en()) return edgez[i+1][j+1][k  ];

  if(dir == wb()) return edgey[i  ][j  ][k  ];
  if(dir == wt()) return edgey[i  ][j  ][k+1];
  if(dir == eb()) return edgey[i+1][j  ][k  ];
  if(dir == et()) return edgey[i+1][j  ][k+1];
        
  if(dir == sb()) return edgex[i  ][j  ][k  ];
  if(dir == st()) return edgex[i  ][j  ][k+1];
  if(dir == nb()) return edgex[i  ][j+1][k  ];
  if(dir == nt()) return edgex[i  ][j+1][k+1];

  if(dir == w()) return faceval[Comp::i()][i  ][j][k];
  if(dir == e()) return faceval[Comp::i()][i+1][j][k];
  if(dir == s()) return faceval[Comp::j()][i][j  ][k];
  if(dir == n()) return faceval[Comp::j()][i][j+1][k];
  if(dir == b()) return faceval[Comp::k()][i][j][k  ];
  if(dir == t()) return faceval[Comp::k()][i][j][k+1];

  return nodeval[i][j][k];  
}

