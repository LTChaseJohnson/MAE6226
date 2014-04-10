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

#Defining panel class
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                       #First endpoint of panel
        self.xb,self.yb = xb,yb                       #Second endpoint of panel
        self.xc,self.yc = sqrt((xb-xa)**2+(yb-ya)**2) #Control center point
        
        #Sets the orientation of the panel
        if (xb-xa<=0.): self.beta = acos((yb-a)/self.length)
        elif (xb-xa>0/): self.beta = pi+acos(-(yb-a)/self.length)
        
        #Sets the location of the panel
        if (self.beta<=pi): self.loc = 'Top Surface'
        else: self.loc = 'Bottom Surface'
        
        self.sigma = 0.                             #Creating initial source
        self.vt = 0.                                #Creating initial tangent velocity
        self.Cp = 0.                                #Creating initial pressure coeff
        



plt.show()