import numpy as np
from  scipy.integrate import quad

deltat = 1.25;
rhov =0.597; 
rhol =958.4;
L =2258e3;
cpl =4215.9;
cpv =2030;
#b0 = np.sqrt(3/np.pi)*(deltat/((rhov/rhol)*L/cpl+(cpl-cpv)/cpl * deltat));
#print(b0);

def integrand(x,b):
  return np.exp(-b**2*((1-x)**(-2)-2*(1-rhov/rhol)*x-1));

def fval(b):
  return rhol*cpl*deltat/(rhov*(L+(cpl-cpv)*deltat)) - 2*b**2 * (quad(integrand,0,1,args=b))[0]

eps = 1e-5
bm = 0;
bp = 10;
fm = fval(bm)
fp = fval(bp)

for i in range(0,100):
#while np.abs(fm-fp) > eps: 
  bc = (bm+bp)/2
  fc = fval(bc)
  print(bc,fc,"|",bm,bp)
  if fc>0:
    bm = bc
    fm = fc
  else:
    bp = bc
    fm = fc

#print(bc)
#print(fc)

