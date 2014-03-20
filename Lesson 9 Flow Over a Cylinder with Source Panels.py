import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

Uinf = input('Enter Freestream velocity: ')

# Defining a cylinder
R = input('Enter cylinder radius: ')
theta = np.linspace(0,(2*pi),100)
xCylinder,yCylinder = R*np.cos(theta),R*np.sin(theta)

# Plotting our cylinder
size = 4
plt.figure(num=1,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize = 16)
plt.ylabel('y',fontsize = 16)
plt.plot(xCylinder,yCylinder,c='b',ls='-',lw=2)
plt.xlim(-R-.1,R+.1)
plt.ylim(-R-.1,R+.1)

# Defining Panel class to calculate geometry for all panels
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya
        self.xb,self.yb = xb,yb
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2
        self.length = np.sqrt((xb-xa)**2+(yb-ya)**2)
        
        # Panel orientation
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi + acos((yb-ya)/self.length)
        
        self.sigma = 0.            # Initializes Sheet Strength
        self.Vt = 0.               # Initializes Tangential Velocity
        self.Cp = 0.               # Initializes Pressure Coefficient

# Building Panels
Np = input('Enter number of panels: ')

xb = R*np.cos(np.linspace(0,(2*pi),Np+1))  # Sets panel endpoints
yb = R*np.sin(np.linspace(0,(2*pi),Np+1))

panel = np.empty(Np,dtype=object)                  # Creating empty panel shell
for i in range(Np):
    panel[i] = Panel(xb[i],yb[i],xb[i+1],yb[i+1])  # Storing panel information
    
# Creating Plot of Panels
size = 6
plt.figure(num=2,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Cylinder Panels')
plt.plot(xb,yb,c='#FF0000',ls='-',lw=2)
plt.plot(xCylinder,yCylinder,c='#00FF00', ls='-',lw=1)
plt.scatter([p.xa for p in panel],[p.ya for p in panel],c='#FF0000',s=40)
plt.scatter([p.xc for p in panel],[p.yc for p in panel],c='#00FF00',s=40)
plt.legend(['Cylinder','Panels','End Points','Center Points'],loc='best',prop={'size':16})
plt.xlim(-R-.1,R+.1)
plt.ylim(-R-.1,R+.1)

# Evaluating integral for flow boundary conditions
def I(pi,pj):
    def func(s):
        return ((pi.xc-(pj.xa-sin(pj.beta)*s))*cos(pi.beta)+(pi.yc-(pj.ya+cos(pj.beta)*s))*sin(pi.beta))\
        /((pi.xc-(pj.xa-sin(pj.beta)*s))**2+(pi.yc-pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]
    
# Setting up System of Equaitons as vectors
A = np.empty((Np,Np),dtype=float)
for i in range(Np):
    for j in range(Np):
        if (i!=j):
            A[i,j] = 0.5/pi * I(panel[i],panel[j])
        else:
            A[i,j] = 0.5

b = =Uinf*np.cos([p.beta for p in panel])

# Solving Vector Equations
var = np.linalg.solve(A,b)
for i in range(len(panel)):
    panel[i].sigma = var[i]

plt.show()