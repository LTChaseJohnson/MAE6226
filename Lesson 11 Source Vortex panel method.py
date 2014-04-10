# Importing Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

# Defining our airfoil from imported geometry
coords = np.loadtxt(fname='C:/Users/akashdhruv/Downloads/s1223.dat')
xp,yp = coords[:,0],coords[:,1]

# Creating our airfoil
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize = 16)
plt.ylabel('y',fontsize = 16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.show()