import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 200                           # Number of points in each direction
xStart,xEnd = -4.0,4.0            # x-direction boundaries
yStart,yEnd = -2.0,2.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid

Uinf = 1.0                        # Free stream velocity
alpha = 0                       # Angle of attack

# Computing the velocity components (u,v)
uFreestream = Uinf*cos(alpha*pi/180)*np.ones((N,N),dtype=float)
vFreestream = Uinf*sin(alpha*pi/180)*np.ones((N,N),dtype=float)

# Plotting the Freestream
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uFreestream,vFreestream,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')

# Merging u and v into Freestream Potential

psiFreestream = Uinf*(Y * cos(alpha*pi/180)- X * sin(alpha*pi/180))

# Creating a function to compute velocity field
def getvelocity(strength,xs,ys,X,Y):
    u = strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v

# Creating a function to compute stream function
def getstreamfunction(strength,xs,ys,X,Y):
    psi = strength/(2*pi) * np.arctan2((Y-ys),(X-xs))
    return psi

# Using functions to calculate one source
sourcestrength1 = 15.0          #Source strength
xSource,ySource = -1.0,0.0     #Source location

uSource,vSource = getvelocity(sourcestrength1,xSource,ySource,X,Y)
psiSource = getstreamfunction(sourcestrength1,xSource,ySource,X,Y)

# Using functions to calculate one sink
sinkstrength1 = -5.0        #Sink strength
xSink,ySink = 1.0,0.0        #Sink location

uSink,vSink = getvelocity(sinkstrength1,xSink,ySink,X,Y)
psiSink = getstreamfunction(sinkstrength1,xSink,ySink,X,Y)

# Plotting the Source
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uSource,vSource,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')

# Plotting the Sink
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uSink,vSink,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSink,ySink,c='#000000',s=80,marker='o')

# Superposition of free stream on source/sink velocity field
u = uFreestream + uSource + uSink
v = vFreestream + vSource + vSink
psi = psiFreestream + psiSource + psiSink

# Plotting combined freestream and source/sink velocity fields
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')
plt.scatter(xSink,ySink,c='#000000',s=80,marker='o')

# Calculating stagnation point
xStagnation = xSource - (sourcestrength1*cos(alpha*pi/180))/(2*pi*Uinf) + xSink - (sinkstrength1*cos(alpha*pi/180))/(2*pi*Uinf)
yStagnation = ySource - (sourcestrength1*sin(alpha*pi/180))/(2*pi*Uinf) + ySink - (sinkstrength1*sin(alpha*pi/180))/(2*pi*Uinf)

if (alpha==0.0):
    plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid') 

# Calculating pressure coefficient
Cp = 1.0-(u**2+v**2)/Uinf**2
print 'Cp = ',Cp

# Pressure strength plot
size = 10
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter([xSource,xSink],[ySource,ySink],c='#CD2305',s=80,marker='o')
plt.contour(X,Y,psi,\
            levels=[0.0],\
            colors='#CD2305',linewidths=2,linestyles='solid')

plt.show()