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

# Defining panel class
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                       #First endpoint of panel
        self.xb,self.yb = xb,yb                       #Second endpoint of panel
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2         #Control center point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)
        
        #Sets the orientation of the panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        #Sets the location of the panel
        if (self.beta<=pi): self.loc = 'Top Surface'
        else: self.loc = 'Bottom Surface'
        
        self.sigma = 0.                             #Creating initial source
        self.vt = 0.                                #Creating initial tangent velocity
        self.Cp = 0.                                #Creating initial pressure coeff
        
# Creating discrete panels for imported geometry
def definePanels(N,xp,yp):
    R = (max(xp)-min(xp))/2
    xCenter = (max(xp)+min(xp))/2
    xCircle = xCenter +R*np.cos(np.linspace(0,2*pi,N+1))
    
    x = np.copy(xCircle)
    y = np.empty_like(x)
    
    xp,yp = np.append(xp,xp[0]),np.append(yp,yp[0])
    
    I = 0
    for i in range(N):
        while (I<len(xp)-1):
            if (xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[I+1]-yp[I])/(xp[I+1]-xp[I])
        b = yp[I+1]-a*xp[I+1]
        y[i] = a*x[i]+b
    y[N] = y[0]
    
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
    
    return panel

N = input('Enter number of panels: ')
panel = definePanels(N,xp,yp)

valX,valY =0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),np.append([p.ya for p in panel],panel[0].ya),'r-',linewidth=1,marker='o',markersize=6)

plt.show()