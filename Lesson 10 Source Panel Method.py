import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

coords = np.loadtxt(fname='C:/Users/chasevjohnson/Downloads/n0012.dat')
xp,yp = coords[:,0],coords[:,1]

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
plt.plot(xp,yp,'k-',linewidth = 2)
plt.axis("equal")

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

valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize = 16)
plt.ylabel('y',fontsize = 16)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa), np.append([p.ya for p in panel],panel[0].ya),linestyle = '-', linewidth = 1, marker = 'o', markersize = 6, color = '#FF0000')
plt.axis("equal")

class Freestream:
	def __init__(self,Uinf,alpha):
		self.Uinf = Uinf                   
		self.alpha = alpha*pi/180

Uinf = input('Enter Freestream velocity: ')
alpha = input('Enter angle of attack: ')
freestream = Freestream(Uinf,alpha)

def I(xci,yci,pj,dxdz,dydz):
	def func(s):
		return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz+(yci-(pj.ya+cos(pj.beta)*s))*dydz)/((xci-(pj.xa-sin(pj.beta)*s))**2 + (yci-(pj.ya+cos(pj.beta)*s))**2)
	return integrate.quad(lambda s:func(s),0.,pj.length)[0]

def buildMatrix(p):
    N = len(p)
    A = np.empty((N,N),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range(N):
        for j in range(N):
            if (i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],cos(p[i].beta),sin(p[i].beta))
    return A

def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N,dtype=float)
    for i in range(N):
        B[i] = -fs.Uinf*cos(fs.alpha-p[i].beta)
    return B

A = buildMatrix(panel)
B = buildRHS(panel,freestream)

var = np.linalg.solve(A,B)
for i in range(len(panel)):
	panel[i].sigma = var[i]

def getTangentVelocity(p,fs):
	N = len(p)
	A = np.zeros((N,N),dtype=float)
	for i in range(N):
		for j in range(N):
			if (i!=j):
				A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),cos(p[i].beta))
	B = fs.Uinf*np.sin([fs.alpha-pp.beta for pp in p])
	var = np.array([pp.sigma for pp in p])
	vt = np.dot(A,var)+B
	for i in range(N):
		p[i].vt = vt[i]	
			
getTangentVelocity(panel,freestream)
     
def getPressureCoefficient(p,fs):
    for i in range(len(p)):
        p[i].Cp = 1-(p[i].vt/fs.Uinf)**2

getPressureCoefficient(panel,freestream)

valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
Cpmin,Cpmax = min([p.Cp for p in panel]),max([p.Cp for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='Top'],[p.Cp for p in panel if p.loc=='Top'],\
		'ro-',linewidth=2)
plt.plot([p.xc for p in panel if p.loc=='Bottom'],\
		[p.Cp for p in panel if p.loc=='Bottom'],\
		'bo-',linewidth=1)
plt.legend(['Top','Bottom'],'best',prop={'size':14})
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.gca().invert_yaxis()
plt.title('Number of panels : %d'%len(panel));
plt.axis("equal")

def getVelocityField(panel,freestream,X,Y):
    Nx,N = X.shape
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.Uinf*cos(freestream.alpha)+0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,1,0) for p in panel])
            v[i,j] = freestream.Uinf*sin(freestream.alpha)+0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,0,1) for p in panel])
    return u,v

Nx,Ny = 80,80
valX,valY = 1.0,2.0
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
X,Y = np.meshgrid(np.linspace(xStart,xEnd,Nx),np.linspace(yStart,yEnd,Ny))

u,v = getVelocityField(panel,freestream,X,Y)

size = 12
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=4,linewidth=1,arrowsize=1,arrowstyle='->')
plt.fill([p.xa for p in panel],[p.ya for p in panel], 'ko-', linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of velocity field')
plt.axis("equal")

Cp = 1.0-(u**2+v**2)/freestream.Uinf**2

size = 12
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
contf = plt.contourf(X,Y,Cp,levels = np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label('$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.fill([p.xa for p in panel],[p.ya for p in panel], 'ko-', linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of Pressure Field')
plt.axis("equal")

plt.show()