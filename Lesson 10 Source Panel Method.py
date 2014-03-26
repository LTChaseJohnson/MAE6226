import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

coords = np.loadtxt('C:/Users/Chase/Dropbox/Graduate_School/AeroDynamics/Python/n0012.dat')
xp,yp = coords[:,0],coords[:,1]

valX,valY = 01,0.2
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
plt.plot(xp,yp,'k-',linewidth = 2)

plt.show()