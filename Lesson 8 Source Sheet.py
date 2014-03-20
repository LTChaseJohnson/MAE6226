import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import integrate

# Creating our mesh grid
N = 200                           # Number of points in each direction
xStart,xEnd = -1.0,1.0            # x axis boundaries
yStart,yEnd = -1.5,1.5            # y axis boundaries
x = np.linspace(xStart,xEnd,N)    # x vector
y = np.linspace(yStart,yEnd,N)    # y vector

X,Y = np.meshgrid(x,y)            # meshgrid of x and y vectors

# Creating Freestream Field
Uinf = input('Enter Freestream velocity: ')

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)

# Creating class to calculate aspects of a Source
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    def velocity(self,X,Y):
        self.u = (self.strength/(2*pi))*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = (self.strength/(2*pi))*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    def streamfunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))

# Defining our FINITE source sheet
Nsources = input('Enter number of sources: ')
strength = input('Enter finite sheet strength: ')
sourcestrength = strength/Nsources
xSource = 0.0                              # Places source sheet on y-axis
ySource = np.linspace(-1.0,1.0,Nsources)   # Evenly spaces sources in sheet

# Creating source sheet
sources = np.empty(Nsources,dtype=object)   # Empty array to fill
for i in range(Nsources):
    sources[i] = Source(sourcestrength,xSource,ySource[i])
    sources[i].velocity(X,Y)

# Combining flow fields using super position
u = uFreestream.copy()
v = vFreestream.copy()                    # I don't know why we are doing this

for s in sources:                        # Adds calculated source fields to
    u = np.add(u,s.u)                    # Existing freestream field
    v = np.add(v,s.v)

# Creating the plot
size = 6
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Finite Source Sheet')
plt.streamplot(X,Y,u,v, density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource*np.ones(Nsources,dtype=float), ySource, c='#FF0000', s=80, marker='o')
velocity = plt.contourf(X,Y,np.sqrt(u**2+v**2), levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(velocity, ticks=[0.0,0.05,0.1],orientation='horizontal')
cbar.set_label(r'Velocity Magnitude $= \sqrt{u^2+v^2}$',fontsize=16)

# Defining our INFINITE source sheet
sigma = input('Enter infinite sheet strength: ')   # Defines strength of sheet
ymin,ymax = -1.0,1.0                               # Defines limits of sheet

uPanel = np.empty((N,N),dtype=float)
vPanel = np.empty((N,N),dtype=float)

# Calculating velocity fields
for i in range(N):
    for j in range(N):
        integrandu = lambda s : X[i,j]/(X[i,j]**2+(Y[i,j]-s)**2)
        uPanel[i,j] = sigma/(2*pi)*integrate.quad(integrandu,ymin,ymax)[0]
        integrandv = lambda s: (Y[i,j]-s)/(X[i,j]**2 + (Y[i,j]-s)**2)
        vPanel[i,j] = sigma/(2*pi)*integrate.quad(integrandv,ymin,ymax)[0]
        
# Combining with freestream using super position principle
u2 = np.add(uFreestream,uPanel)
v2 = np.add(vFreestream,vPanel)

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Infite Source Sheet')
plt.streamplot(X,Y,u2,v2, density=2, linewidth=1, arrowsize=1, arrowstyle='->')
plt.axvline(0.0, (ymin-yStart)/(yEnd-yStart), (ymax-yStart)/(yEnd-yStart), color='#FF0000', linewidth=4)
velocity = plt.contourf(X,Y,np.sqrt(u2**2+v2**2), levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(velocity, ticks=[0.0,0.05,.01], orientation='horizontal')
cbar.set_label(r'Velocity Magnitude $= \sqrt{u^2+v^2}$')
plt.show()
