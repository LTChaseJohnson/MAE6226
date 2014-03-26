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

class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya
        self.xb,self.yb = xb,yb
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)
        
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        if (self.beta<=pi): self.loc = 'Top'
        else: self.loc = 'Bottom'
        
        self.sigma = 0.
        self.vt = 0.
        self.Cp = 0.

def definePanels(N,xp,yp):
    R = (max(xp)-min(xp))/2
    xc,yc = (max(xp)+min(xp))/2,(max(yp)+min(yp))/2
    xCircle = xc + R*np.cos(np.linspace(0,2*pi,N+1))
    yCircle = yc + R*np.sin(np.linspace(0,2*pi,N+1))
    
    x = np.copy(xCircle[0:-1])
    y = np.empty_like(x)
    
    I = 0
    for i in range(N):
        while (I<len(xp)-1):
            if (xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[(I+1)%len(yp)]-yp[I])/(xp[(I+1)%len(yp)]-xp[I])
        b = yp[(I+1)%len(yp)]-a*xp[(I+1)%len(xp)]
        y[i] = a*x[i]+b
        
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[(i+1)%N],y[(i+1)%N])
        
    return panel

N = input('Enter number of panels: ')
panel = definePanels(N,xp,yp)

xmin2,xmax2 = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin2,ymax2 = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart2,xEnd2 = xmin2-valX*(xmax2-xmin2),xmax2+valX*(xmax2-xmin2)
yStart2,yEnd2 = ymin2-valY*(ymax2-ymin2),ymax2+valY*(ymax2-ymin2)
plt.figure(figsize=(size,(yEnd2-yStart2)/(xEnd2-xStart2)*size))
plt.grid(True)
plt.xlabel('x',fontsize = 16)
plt.ylabel('y',fontsize = 16)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa), np.append([p.ya for p in panel],panel[0].ya),linestyle = '-', linewidth = 1, marker = 'o', markersize = 6, color = '#FF0000')

plt.show()