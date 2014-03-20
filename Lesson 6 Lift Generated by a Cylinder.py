import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 200                                                   # Number of points in each direction
xStart,xEnd = -2.0,2.0                                    # Sets the limits of x vector
yStart,yEnd = -1.0,1.0                                    # and y vector
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

kappa = input('Doublet strength: ')
xDoublet,yDoublet = input('Doublet location x,y: ')
# Developing external free stream based on user input
Uinf = input('Enter freestream velocity: ')
alpha1 = input('Enter freestream angle of attack: ')
alpha = alpha1*(pi/180)
uFreestream = Uinf*np.cos(alpha)*np.ones((N,N),dtype=float) #Creates usable array
vFreestream = Uinf*np.sin(alpha)*np.ones((N,N),dtype=float) #Creates usable array
psiFreestream = Uinf*np.cos(alpha)*Y

def getfunctionsdoublet(strength,xd,yd,X,Y):
    u = -strength/(2*pi)*((xd-X)**2-(yd-Y)**2)/(((xd-X)**2+(yd-Y)**2)**2)
    v = -strength/(2*pi)*(2*(xd-X)*(yd-Y))/(((xd-X)**2+(yd-Y)**2)**2)
    psi = -strength/(2*pi)*(yd-Y)/((xd-X)**2+(yd-Y)**2)
    return u,v,psi
    
uDoublet,vDoublet,psiDoublet = getfunctionsdoublet(kappa,xDoublet,yDoublet,X,Y)

u1 = uDoublet + uFreestream
v1 = vDoublet + vFreestream
psi1 = psiDoublet + psiFreestream

# Determining stagnation points for general flow
xStag1,yStag2 =-sqrt(kappa/(2*pi*Uinf))+xDoublet,sin(alpha)*sqrt(kappa/(2*pi*Uinf))+yDoublet
xStag2,yStag1 =sqrt(kappa/(2*pi*Uinf))+xDoublet,sin(alpha+pi)*sqrt(kappa/(2*pi*Uinf))+yDoublet

# Determining radius
R = sqrt(kappa/(2*pi*Uinf))

# Plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Cylinder in a Freestream')
plt.ylim(yStart,yEnd)
plt.xlim(xStart,xEnd)
plt.streamplot(X,Y,u1,v1,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
circle = plt.Circle((xDoublet,yDoublet),radius=R,color='#ff0000')
plt.gca().add_patch(circle)
plt.scatter(xDoublet,yDoublet,c='#CDFFFF')
plt.scatter([xStag1,xStag2],[yStag1,yStag2],c='#400D12',s=80,marker='o')

gammaprime = 4*pi*Uinf*R
print ('4*pi*Uinf*R= '),gammaprime
gamma = input('Vortex strength: ')
xVortex,yVortex = xDoublet,yDoublet

def getfunctionsvortex(strength,xv,yv,X,Y):
    u = strength/(2*pi)*((yv-Y)/((xv-X)**2+(yv-Y)**2))
    v = -strength/(2*pi)*((xv-X)/((xv-X)**2+(yv-Y)**2))
    psi = strength/(4*pi)*np.log((xv-X)**2+(yv-Y)**2)
    return u,v,psi

uVortex,vVortex,psiVortex = getfunctionsvortex(gamma,xVortex,yVortex,X,Y)

u = uVortex + uDoublet + uFreestream
v = vVortex + vDoublet + vFreestream

# New Stagnation Points
rs = sqrt(kappa/(2*pi*Uinf*cos(alpha)))
thetas = (-gamma/(2*pi*rs))/((1/(Uinf*cos(alpha)))+(Uinf*cos(alpha)))
xstag1,ystag1 = rs*np.cos(thetas),-rs*np.sin(thetas)
xstag2,ystag2 = -rs*np.cos(thetas),-rs*np.sin(thetas)


plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Cylinder in a Freestream with Vortex')
plt.ylim(yStart,yEnd)
plt.xlim(xStart,xEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
circle = plt.Circle((xDoublet,yDoublet),radius=R,color='#ff0000')
plt.gca().add_patch(circle)
plt.scatter(xDoublet,yDoublet,c='#CDFFFF')
plt.scatter([xstag1,xstag2],[ystag1,ystag2],c='#400D12',s=80,marker='o')

theta = np.linspace(0,2*pi,200)
C1 = np.cos(alpha)*np.cos(theta)+np.sin(alpha)*np.sin(theta)
C2 = np.sin(alpha)*np.cos(theta)-np.cos(alpha)*np.sin(theta)
utheta = Uinf*(C1+C2)+gamma/(2*pi*R)
uthetan = Uinf*(C1+C2)
Cp = 1-(utheta/Uinf)**2
Cpn = 1-(uthetan/Uinf)**2

plt.figure(figsize=(size,size))
plt.grid(True)
plt.xlabel(r'$\theta$',fontsize=18)
plt.ylabel(r'$C_p$',fontsize=18)
plt.xlim(theta.min(),theta.max())
plt.plot(theta,Cp,color='#CD2305',linewidth=2,linestyle='-')
plt.plot(theta,Cpn,color='g',linewidth=2,linestyle='-')
plt.legend(['with vortex','without vortex'],loc='best',prop={'size':16});

# Plotting pressure contour
Cp1 = 1.0-(u1**2+v1**2)/Uinf**2
Cpv = 1.0-(u**2+v**2)/Uinf**2

plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Pressure Gradient with no Vortex')
contf=plt.contourf(X,Y,Cp1,levels=np.linspace(-5.0,2.0,100),extend='both')
cbar=plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')

plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Pressure Gradient with Vortex')
contf=plt.contourf(X,Y,Cpv,levels=np.linspace(-5.0,2.0,100),extend='both')
cbar=plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')

plt.show()