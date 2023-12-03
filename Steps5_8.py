# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:48:26 2023

@author: emrul
"""

import numpy as np
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

def TwoD_LinearConvection(nx,ny,speed):
    
    # Assign Variables
    dx = 2/(nx-1)
    dy = 2/(ny-1)
    nt = 100
    c = speed
    sigma = 0.2
    dt = sigma*dx/c
    
    x = np.linspace(0,2,nx)
    y = np.linspace(0,2,ny)
    
    u = np.ones((ny,nx))
    un = np.ones((ny,nx))
    
    # Assign Initial Conditions
    
    u[int(0.5/dy):int(1/dy+1),int(.5/dx):int(1/dx+1)] = 2
    
    print(u)
    
    #Plot Initial Conditions
    """
    fig = pyplot.figure(figsize=(11,7),dpi=100)
    ax = fig.add_subplot(projection = '3d')
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,u[:],cmap=cm.viridis)
    pyplot.show()
    """
    """
    for n in range(nt+1):
        un = u.copy()
        row,col = u.shape

        for j in range(1,row):
            for i in range(1,col):

                u[j,i] = un[j,i] - (c*(dt/dx)*(un[j,i]-un[j,i-1])) - (c*(dt/dy)*(un[j,i]-un[j-1,i]))

                u[0,:]  = 1
                u[-1,:] = 1
                u[:,0]  = 1
                u[:,-1] = 1
    
    """
    
    for n in range(nt+1):
        un = u.copy()

        u[1:,1:] = un[1:,1:] - (c*(dt/dx)*(un[1:,1:]-un[1:,:-1])) - (c*(dt/dy)*(un[1:,1:]-un[:-1,1:]))

        u[0,:]  = 1
        u[-1,:] = 1
        u[:,0]  = 1
        u[:,-1] = 1

    fig = pyplot.figure(figsize=(11,7),dpi=100)
    ax = fig.add_subplot(projection = '3d')
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,u[:],cmap=cm.viridis)
    pyplot.show()

def TwoD_Convection(ny,nx):

    # Assign Variables

    dx = 2/(nx-1)
    dy = 2/(ny-1)
    nt = 100
    sigma = 0.2
    dt = sigma*dx
    
    x = np.linspace(0,2,nx)
    y = np.linspace(0,2,ny)
    
    u = np.ones((ny,nx))
    un = np.ones((ny,nx))
    v = np.ones((ny,nx))
    vn = np.ones((ny,nx))

    # Boundary & Initial Conditions
    
    u[int(0.5/dy):int(1/dy+1),int(.5/dx):int(1/dx+1)] = 2
    v[int(0.5/dy):int(1/dy+1),int(.5/dx):int(1/dx+1)] = 2

    # Iterate Over Time

    for n in range(nt+1):

        un = u.copy()
        vn = v.copy()

        u[1:,1:] = un[1:,1:] - (un[1:,1:]*(dt/dx)*(un[1:,1:]-un[1:,:-1])) - (vn[1:,1:]*(dt/dy)*(un[1:,1:]-un[:-1,1:]))
        v[1:,1:] = vn[1:,1:] - (un[1:,1:]*(dt/dx)*(vn[1:,1:]-vn[1:,:-1])) - (vn[1:,1:]*(dt/dy)*(vn[1:,1:]-vn[:-1,1:]))

        u[0,:]  = 1
        u[-1,:] = 1
        u[:,0]  = 1
        u[:,-1] = 1

        v[0,:]  = 1
        v[-1,:] = 1
        v[:,0]  = 1
        v[:,-1] = 1

    fig = pyplot.figure(figsize=(11,7),dpi=100)
    ax = fig.add_subplot(projection = '3d')
    X,Y = np.meshgrid(x,y)
    ax.plot_surface(X,Y,u,cmap=cm.viridis, rstride=2,cstride=2)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

    fig2 = pyplot.figure(figsize=(11,7),dpi=100)
    ax2 = fig2.add_subplot(projection = '3d')
    ax2.plot_surface(X,Y,v,cmap=cm.viridis, rstride=2,cstride=2)
    ax2.set_xlabel('$x$')
    ax2.set_ylabel('$y$')

    pyplot.show()

def TwoD_Diffusion(nx,ny,nt):

    # Assign Variables
    dx = 2/(nx-1)
    dy = 2/(ny-1)
    nu=.05
    sigma = 0.25
    dt = sigma*dx*dy/nu
    
    x = np.linspace(0,2,nx)
    y = np.linspace(0,2,ny)
    
    u = np.ones((ny,nx))
    un = np.ones((ny,nx))
    
    # Assign Initial Conditions
    
    u[int(0.5/dy):int(1/dy+1),int(.5/dx):int(1/dx+1)] = 2

    # Iterate Over Time

    for n in range(nt+1):
        un = u.copy()

        u[1:-1,1:-1] = un[1:-1,1:-1] + (nu*(dt/dx**2)*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])) + (nu*(dt/dy**2)*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]))

        u[0,:]  = 1
        u[-1,:] = 1
        u[:,0]  = 1
        u[:,-1] = 1

    fig = pyplot.figure(figsize=(11,7),dpi=100)
    ax = fig.add_subplot(projection = '3d')
    X,Y = np.meshgrid(x,y)
    ax.plot_surface(X,Y,u[:],cmap=cm.viridis)
    ax.set_zlim(1, 2.5)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    pyplot.show()

def TwoD_BurgersEqu(nx,ny):

    # Define Variables

    dx = 2/(nx-1)
    dy = 2/(ny-1)
    nt=120
    nu=.001
    sigma = 0.009
    dt = sigma*dx*dy/nu
    
    x = np.linspace(0,2,nx)
    y = np.linspace(0,2,ny)
    
    u = np.ones((ny,nx))
    un = np.ones((ny,nx))
    v = np.ones((ny,nx))
    vn = np.ones((ny,nx))
    
    # Define BC's and IC's

    u[int(0.5/dy):int(1/dy+1),int(.5/dx):int(1/dx+1)] = 2
    v[int(0.5/dy):int(1/dy+1),int(.5/dx):int(1/dx+1)] = 2

    #Iterate Over Time

    for n in range(nt+1):
        un = u.copy()
        vn = v.copy()
        # Equation for variation of x component

        u[1:-1,1:-1] = un[1:-1,1:-1] - (un[1:-1,1:-1]*(dt/dx)*(un[1:-1,1:-1]-un[1:-1,0:-2])) - (vn[1:-1,1:-1]*(dt/dy)*(un[1:-1,1:-1]-un[0:-2,1:-1])) + (nu*(dt/(dx**2))*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])) + (nu*(dt/(dy**2))*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]))

        # Equation for variation of y component

        v[1:-1,1:-1] = vn[1:-1,1:-1] - (un[1:-1,1:-1]*(dt/dx)*(vn[1:-1,1:-1]-vn[1:-1,0:-2])) - (vn[1:-1,1:-1]*(dt/dy)*(vn[1:-1,1:-1]-vn[0:-2,1:-1])) + (nu*(dt/(dx**2))*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])) + (nu*(dt/(dy**2))*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1]))

        u[0,:]  = 1
        u[-1,:] = 1
        u[:,0]  = 1
        u[:,-1] = 1

        v[0,:]  = 1
        v[-1,:] = 1
        v[:,0]  = 1
        v[:,-1] = 1
    
    # Plot the Results

    fig = pyplot.figure(figsize=(11,7),dpi=100)
    ax = fig.add_subplot(projection = '3d')
    X,Y = np.meshgrid(x,y)
    ax.plot_surface(X,Y,u,cmap=cm.viridis, rstride=2,cstride=2)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

    fig2 = pyplot.figure(figsize=(11,7),dpi=100)
    ax2 = fig2.add_subplot(projection = '3d')
    ax2.plot_surface(X,Y,v,cmap=cm.viridis, rstride=2,cstride=2)
    ax2.set_xlabel('$x$')
    ax2.set_ylabel('$y$')

    pyplot.show()

if __name__ == '__main__':
    
    #TwoD_LinearConvection(101,101,2)
    #TwoD_Convection(101,101)
    #TwoD_Diffusion(31,31,50)
    TwoD_BurgersEqu(81,81)