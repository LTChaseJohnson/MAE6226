import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50                                               # Number of points in each direction
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)                       # Creates linear array for x-axis
y = np.linspace(yStart,yEnd,N)                       # Creates linear array for y-axis

X,Y = np.meshgrid(x,y)                               # Creates a mesh grid using x and y

k = input('Doublet Strength: ')                      # Doublet strength
xDoublet,yDoublet = input('Doublet Location: x,y ')  # Doublet location

# Creates function getvelocity that determines velocity field based
# on user input of strength and doublet location
def getvelocity(strength,xd,yd,X,Y):
    u = -strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = -strength/(2*pi)*(2*(X-xd)*(Y-yd))/((X-xd)**2+(Y-yd)**2)**2
    return u,v

# Creates function getstreamfunction that determines potential field
# based on user input of strength and doublet location
def getstreamfunction(strength,xd,yd,X,Y):
    psi = -strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi

# Calculating velocity and stream function
uDoublet,vDoublet = getvelocity(k,xDoublet,yDoublet,X,Y)
psiDoublet = getstreamfunction(k,xDoublet,yDoublet,X,Y)

# Plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size)) #Scales the figure based on user input
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Velocity field of a Doublet')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uDoublet,vDoublet,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')  #Plots the doublet as large red dot

# Developing external free stream based on user input
Uinf = input('Enter freestream velocity: ')
alpha1 = input('Enter freestream angle of attack: ')
alpha = alpha1*(pi/180)
uFreestream = Uinf*cos(alpha)*np.ones((N,N),dtype=float) #Creates usable array
vFreestream = Uinf*sin(alpha)*np.ones((N,N),dtype=float) #Creates usable array
psiFreestream = Uinf*Y

# Combining Freestream with Doublet Flow
u = uFreestream + uDoublet
v = vFreestream + vDoublet
psi = psiFreestream + psiDoublet

# Determining stagnation points for general flow
xStag1,yStag2 =-sqrt(k/(2*pi*Uinf))+xDoublet,sin(alpha)*sqrt(k/(2*pi*Uinf))+yDoublet
xStag2,yStag1 =sqrt(k/(2*pi*Uinf))+xDoublet,sin(alpha+pi)*sqrt(k/(2*pi*Uinf))+yDoublet

# Plotting combined flow
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Combined Freestream and Doublet Flow')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1.0,arrowsize=1,arrowstyle='->')
plt.contour(X,Y,psi,levels=[0.0],colors='#cd2305',linewidth=2,linestyles='solid')
plt.scatter(xDoublet,yDoublet,c='b',s=80,marker='o')
plt.scatter([xStag1,xStag2],[yStag1,yStag2],c='g',s=80,marker='o')

# Computing pressure coefficient
Cp = 1.0-(u**2+v**2)/Uinf**2

# Plotting pressure contour
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Pressure Gradient')
contf=plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar=plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')
#plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid')
plt.scatter([xStag1,xStag2],[yStag1,yStag2],c='g',s=80,marker='o')
plt.show()