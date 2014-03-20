import numpy as np
import matplotlib.pyplot as plt
from math import *

# Making meshgrid
N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

#Vortex user inputs
gamma = input('+Vortex strength: ')
xVortex,yVortex = input('+Vortex location x,y: ')
ngamma = input('-Vortex strength: ')
xnVortex,ynVortex = input('-Vortex location x,y: ')
Uinf = input('Freestream velocity: ')
alpha1 = input('Freestream angle of attack: ')
alpha = alpha1*(pi/180)

#Freestream contributions
uFreestream = Uinf*cos(alpha)*np.ones((N,N),dtype=float) #Creates usable array
vFreestream = Uinf*sin(alpha)*np.ones((N,N),dtype=float) #Creates usable array
psiFreestream = Uinf*Y

#Velocity function
def getvelocityvortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

#Stream function
def getstreamfunctionvortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi
    
#Determining fields
uVortex,vVortex = getvelocityvortex(gamma,xVortex,yVortex,X,Y)
unVortex,vnVortex = getvelocityvortex(ngamma,xnVortex,ynVortex,X,Y)
psiVortex = getstreamfunctionvortex(gamma,xVortex,yVortex,X,Y)
psinVortex = getstreamfunctionvortex(ngamma,xnVortex,ynVortex,X,Y)

u = uVortex + uFreestream + unVortex
v = vVortex + vFreestream + vnVortex
psi = psiVortex + psiFreestream + psinVortex

#Plotting Vortex
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Single Vortex')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uVortex,vVortex,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='#CD2305',s=80,marker='o')

#Plotting Vortex Pair
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Vortex Pair')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uVortex+unVortex,vVortex+vnVortex,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='#CD2305',s=80,marker='o')
plt.scatter(xnVortex,ynVortex,c='#000000',s=80,marker='o')

#Plotting Vortex in Uniform Flow
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Vortex in Uniform Flow')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uVortex+uFreestream,vVortex+vFreestream,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='#CD2305',s=80,marker='o')

#Plotting combined field
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Vortex Pair in Uniform Field')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid')
plt.scatter(xVortex,yVortex,c='#CD2305',s=80,marker='o')
plt.scatter(xnVortex,ynVortex,c='#000000',s=80,marker='o')

# Computing pressure coefficient
Cp = 1.0-(u**2+v**2)/Uinf**2

# Plotting pressure contour
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Pressure Gradient')
contf=plt.contourf(X,Y,Cp,levels=np.linspace(-10.0,2.0,100),extend='both')
cbar=plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-10.0,-5.0,-3.0,-2.0,-1.0,0.0,1.0,2.0])
plt.scatter([xVortex,xnVortex],[yVortex,ynVortex],c='#CD2305',s=80,marker='o')
plt.show()