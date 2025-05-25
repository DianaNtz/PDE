"""Poisson equation solver via spectral methods"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#create grid
def Cheb_grid_and_D(xmin,xmax,Nx):
     x=np.zeros(Nx)
     cx=(xmax-xmin)/2
     dx=(xmax+xmin)/2
     for j in range(0,Nx):
         x[j]=-np.cos(j*np.pi/(Nx-1))*cx+dx
     Dx=np.zeros((Nx,Nx))
     Dx[0][0]=-(2*(Nx-1)**2+1)/6
     Dx[Nx-1][Nx-1]=(2*(Nx-1)**2+1)/6
     for i in range(0,Nx):
         if(i==0 or i==Nx-1):
             ci=2
         else:
             ci=1
         for j in range(0,Nx):
             if(j==0 or j==Nx-1):
                 cj=2
             else:
                 cj=1
             if(i==j and i!=0 and i!=Nx-1):
                Dx[i][i]=-((x[i]-dx)/cx)/(1-((x[i]-dx)/cx)**2)*0.5
             if(i!=j):
                Dx[i][j]=ci/cj*(-1)**(i+j)/(x[i]/cx-x[j]/cx)

     Dx=Dx/cx
     return x,Dx,np.dot(Dx,Dx)
#Initial values
Nx=30
Ny=40

xmin=0
xmax=1.0

ymax=1.0
ymin=0.0


y,Dy,D2y=Cheb_grid_and_D(ymin,ymax,Ny)
x,Dx,D2x=Cheb_grid_and_D(xmin,xmax,Nx)

d2x=np.kron(np.identity((Ny)),D2x)
d2y=np.kron(D2y,np.identity((Nx)))

f=np.zeros(Ny*Nx)


for j in range(0,Ny):
          for i in range(0,Nx):
                  f[i+j*Nx]=-2*np.pi**2*np.sin(np.pi*y[j])*np.sin(np.pi*x[i])
#create Laplace operator
L=d2y+d2x
#setting boundaries
bound=np.kron(np.identity((Ny)),np.identity((Nx)))
b=f
for j in range(0,Ny):
        for i in range(0,Nx):
            if(i==0):
               b[i+j*Nx]=0
            if(i==Nx-1):
               b[i+j*Nx]=0
            if(j==Ny-1):
               b[i+j*Nx]=0
            if(j==0):
               b[i+j*Nx]=0



for k in range(0,Ny*Nx):
        if(k%Nx==0):
           #z bei 0
           L[k]=bound[k]
        if(k=Ny*Nx-Nx):
           L[k]=bound[k]
           #rho bei Nrho
