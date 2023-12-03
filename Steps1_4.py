# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 10:11:18 2023

@author: emrul
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy
from sympy import init_printing
init_printing(use_latex=True)
from sympy.utilities.lambdify import lambdify

def oneD_linconvection(nx,ny,speed):
    
    dx = 2/(nx-1)
    dy = 2/(nx-1)
    nt = 20
    c = speed
    
    sigma = 0.5
    
    dt = sigma*dx/c
    
    u = np.ones(nx)
    
    u[int(0.5/dx):int(1/dx+1)] = 2
    
    un = np.ones(nx)
    
    for n in range(nt):
         un = u.copy()
         
         for i in range(1,nx):
             
             u[i] = un[i] - c*(dt/dx)*(un[i] - un[i-1])
             
    plt.plot(np.linspace(0,2,nx),u)
    plt.show()

def oneD_convection(nx):
    
    dx = 2/(nx-1)
    nt = 130
    
    sigma = 0.5
    
    u = np.ones(nx)
    
    u[int(0.5/dx):int(1/dx+1)] = 2
    
    dt = sigma*dx/max(u)
    
    un = np.ones(nx)
    
    for n in range(nt):
         un = u.copy()
         
         for i in range(1,nx):
             
             u[i] = un[i] - un[i]*(dt/dx)*(un[i] - un[i-1])
             
    plt.plot(np.linspace(0,2,nx),u)
    plt.show()
    
def oneD_diffusion (nx):
    
    dx = 2/(nx-1)
    nt = 20
    nu = 0.3
    sigma = 0.2
    
    dt = sigma*(dx**2)/nu
    
    u = np.ones(nx)
    
    u[int(0.5/dx):int(1/dx+1)] = 2
    
    un = np.ones(nx)
    
    for n in range(nt):
         un = u.copy()
         
         for i in range(1,nx-1):
             
             u[i] = un[i] + nu*(dt/dx**2)*(un[i+1] -2*un[i] + un[i-1])
             
    plt.plot(np.linspace(0,2,nx),u)
    plt.show()
    
def burgers_eqs(nx):
    x, nu, t =sympy.symbols('x nu t') # defining sybolic variables
    
    phi = (sympy.exp(-(x-4*t)**2/(4*nu*(t+1)))+sympy.exp(-(x-4*t-2*sympy.pi)**2/(4*nu*(t+1)))) # writing down the equations with sybollic varibles
    #print(phi)
    phiprime = phi.diff(x) # taking derivative with respect to x
    #print(phiprime)
    u = -2 * nu * (phiprime / phi) + 4 #defining u with symbolic variables
    ufunc = lambdify((t, x, nu), u) # lambdify solves the symbolic equation with given parameters for symbolic variables
    # end of sympy equation construction
    
    
    dx = 2*np.pi/(nx-1)
    nu = .07
    dt = dx*nu
    nt = 2
    print(nt)
    x = np.linspace(0,2*np.pi,nx)
    
    un = np.empty(101)
    
    t = 0
    
    u = np.asarray([ufunc(t,x0,nu)for x0 in x])
    #print(u)
    
    for n in range(nt):
        un = u.copy()
        
        for i in range(1,nx-1):
            u[i] = un[i]-un[i]*(dt/dx)*(un[i]-un[i-1]) + nu*(dt/(dx**2))*(un[i+1]-2*un[i]+un[i-1])
            
        u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 *(un[1] - 2 * un[0] + un[-2])
        u[-1] = u[0]
    
    u_analytical = np.asarray([ufunc(nt * dt, xi, nu) for xi in x])
    
    plt.figure(figsize=(11, 7), dpi=100)
    plt.plot(x,u, marker='o', lw=2, label='Computational')
    plt.plot(x, u_analytical, label='Analytical')
    plt.xlim([0, 2 * np.pi])
    plt.ylim([0, 10])
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':
    
    #oneD_linconvection(101, 1)
    #oneD_convection(101)
    #oneD_diffusion(81)
    burgers_eqs(100)
    
    