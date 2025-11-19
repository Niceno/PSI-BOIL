import numpy as np
import matplotlib.pyplot as plt


blending_angle = 40./180.*np.pi
n0sq = 1.-1./(1.+np.tan(blending_angle)*np.tan(blending_angle))

def bfactor(x):
  return (1/(1+(np.tan(x/180*np.pi))**2)-n0sq)/(1-2*n0sq)

print(bfactor(40),bfactor(40.9),bfactor(41))

x = np.linspace(40,50)
y = bfactor(x)

plt.plot(x,y)
plt.show()
