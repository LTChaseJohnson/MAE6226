import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50                              #Number of points in array

xStart,xEnd = -2.0,2.0              #x-direction boundaries
yStart,yEnd = -1.0,1.0              #y-direction boundaries
 
x = np.linspace(xStart,xEnd,N)      #x linear array with N points
y = np.linspace(yStart,yEnd,N)      #y linear array with N points

print 'x= ',x
print 'y= ',y

X,Y = np.meshgrid(x,y)              #creates X Y as square arrays of N 
                                    #rows of x or y linear arrays
size = 10

# plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))  #sets figure width,height
# plt.xlabel('x',fontsize=16)
# plt.ylabel('y',fontsize=16)
# plt.xlim(xStart,xEnd)
# plt.ylim(yStart,yEnd)

# plt.scatter(X,Y,s=10,c='#CD2305',marker='o',linewidth=0.1)

# Creates a line plot of a source

strengthsource = 5.0               #Sets the source strength

xSource,ySource = -1.0,0.0         #Sets location of the source

uSource = np.empty((N,N),dtype=float)  #Creates empty 2-D array for u
vSource = np.empty((N,N),dtype=float)  #Creates empty 2-D array for v

for i in range(N):
    for j in range(N):
        uSource[i,j] = strengthsource/(2*pi)*(X[i,j]-xSource)/((X[i,j]-xSource)**2+(Y[i,j]-ySource)**2)
        vSource[i,j] = strengthsource/(2*pi)*(Y[i,j]-ySource)/((X[i,j]-xSource)**2+(Y[i,j]-ySource)**2)

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)

plt.streamplot(X,Y,uSource,vSource,density=5.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')

# Creates a line plot of a sink

strengthsink = -5.0               #Sets the sink strength

xSink,ySink = 1.0,0.0             #Sets location of the sink

uSink = np.empty((N,N),dtype=float)  #Creates empty 2-D array for u
vSink = np.empty((N,N),dtype=float)  #Creates empty 2-D array for v

for i in range(N):
    for j in range(N):
        uSink[i,j] = strengthsink/(2*pi)*(X[i,j]-xSink)/((X[i,j]-xSink)**2+(Y[i,j]-ySink)**2)
        vSink[i,j] = strengthsink/(2*pi)*(Y[i,j]-ySink)/((X[i,j]-xSink)**2+(Y[i,j]-ySink)**2)

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)

plt.streamplot(X,Y,uSink,vSink,density=5.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xSink,ySink,c='#000000',s=80,marker='o')

# Creates a plot of source/sink pair

uPair = np.empty_like(uSource)      #Creates an empty array of uSource dimensions
vPair = np.empty_like(vSource)

uPair = uSource + uSink
vPair = vSource + vSink

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uPair,vPair,density=5.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')
plt.scatter(xSink,ySink,c='#000000',s=80,marker='o')
plt.show()