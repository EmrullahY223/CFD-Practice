from matplotlib import pyplot, cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def plot2DLaplace(x,y,p):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    ax = fig.add_subplot(projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    pyplot.show()

def laplace2d(p,y,dx,dy,l1norm_target):


    l1norm = 1
    pn = np.empty_like(p)

    while l1norm > l1norm_target:
        pn = p.copy()

        p[1:-1,1:-1] = ((dy**2)*(pn[1:-1,2:]+pn[1:-1,0:-2]) + (dx**2)*(pn[2:,1:-1]+pn[0:-2,1:-1]))/(2*((dx**2)+(dy**2)))

        p[:,0] = 0  #P=0 @ x = 0
        p[:,-1] = y #P=y @ x = 2 
        p[0,:] = p[1,:]     #dp/dy = 0  @ y=0
        p[-1,:] = p[-2,:]   #dp/dy = 0  @ y=1

        l1norm = (np.sum(np.abs(p[:]) - np.abs(pn[:])) /np.sum(np.abs(pn[:])))
    
    return p

def Poisson2D(p,b,dx,dy,nt):
    l1norm = 1
    pn = np.empty_like(p)

    for it in range(nt):

        pn = p.copy()

        p[1:-1,1:-1] = ((dy**2)*(pn[1:-1,2:]+pn[1:-1,0:-2]) + (dx**2)*(pn[2:,1:-1]+pn[0:-2,1:-1])-(b[1:-1,1:-1]*dx**2*dy**2))/(2*((dx**2)+(dy**2)))

        p[:,0] = 0  #P=0 @ x = 0
        p[:,-1] = 0 #P=y @ x = 2 
        p[0,:] = 0     #dp/dy = 0  @ y=0
        p[-1,:] = 0   #dp/dy = 0  @ y=1

    
    return p

def calc_b(b,rho,dt,u,v,dx,dy):

    b[1:-1, 1:-1] = rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + 
                                     (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) - 
                                     ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 
                                     2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - 
                                     ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)
    #b[1:-1, 1:-1] = rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) - ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *(v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)
    return b

def N_S_Poisson(p,dx,dy,b,nit):

    pn = np.empty_like(p)

    for it in range(50):

        pn = p.copy()

        p[1:-1,1:-1] = ((dy**2)*(pn[1:-1,2:]+pn[1:-1,0:-2]) + 
                        (dx**2)*(pn[2:,1:-1]+pn[0:-2,1:-1])-
                        (b[1:-1,1:-1]*(dx**2)*(dy**2)))/(2*((dx**2)+(dy**2)))

        p[:,-1] = p[:,-2]
        p[:,0] = p[:,1]
        p[0,:] = p[1,:]
        p[-1,:] = 0

        return p
    
def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu,nx,ny ,nit):
    # Variables
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))
    
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        
        b = calc_b(b, rho, dt, u, v, dx, dy)
        p = N_S_Poisson(p, dx, dy, b, nit)
        
        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * dt / dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

        u[0, :]  = 0
        u[:, 0]  = 0
        u[:, -1] = 0
        u[-1, :] = 1    # set velocity on cavity lid equal to 1
        v[0, :]  = 0
        v[-1, :] = 0
        v[:, 0]  = 0
        v[:, -1] = 0
        
        
    return u, v, p

def calc_b_periodic(b,rho,dt,u,v,dx,dy):

    b[1:-1, 1:-1] = rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + 
                                     (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) - 
                                     ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 
                                     2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - 
                                     ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)
    
    # Periodic BC Calculation

    # Periodic @ x = 0

    b[1:-1,0] = rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) + 
                                     (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) - 
                                     ((u[1:-1, 1] - u[1:-1, 0]) / (2 * dx))**2 - 
                                     2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) * (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx)) - 
                                     ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2)
    
    #periodic @ x = 2

    b[1:-1,-1] = rho * (1 / dt * ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx) + 
                                     (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) - 
                                     ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 - 
                                     2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) * (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) - 
                                     ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2)

    return b

def poisson_periodic(p,dx,dy,b,nit):

    pn = np.empty_like(p)

    for it in range(50):

        pn = p.copy()

        p[1:-1,1:-1] = ((dy**2)*(pn[1:-1,2:]+pn[1:-1,0:-2]) + 
                        (dx**2)*(pn[2:,1:-1]+pn[0:-2,1:-1])-
                        (b[1:-1,1:-1]*(dx**2)*(dy**2)))/(2*((dx**2)+(dy**2)))
        
        # Calculate Periodic BC

        # Periodic at x = 0

        p[1:-1,0] = ((dy**2)*(pn[1:-1,1]+pn[1:-1,-1]) + 
                        (dx**2)*(pn[2:,0]+pn[0:-2,0])-
                        (b[1:-1,0]*(dx**2)*(dy**2)))/(2*((dx**2)+(dy**2)))
        
        # Periodic at x = 2
        
        p[1:-1,-1] = ((dy**2)*(pn[1:-1,0]+pn[1:-1,-2]) + 
                        (dx**2)*(pn[2:,-1]+pn[0:-2,-1])-
                        (b[1:-1,-1]*(dx**2)*(dy**2)))/(2*((dx**2)+(dy**2)))
        
        p[0,:] = p[1,:]
        p[-1,:] = p[-2,:]

        return p
    
def channel_flow_periodic(nt, u, v, dt, dx, dy, p, rho, nu,nx,ny ,nit,F):

    # Variables
    udiff = 1
    stepcount = 0
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))
    
    while udiff > 1e-3:
        un = u.copy()
        vn = v.copy()
        
        b = calc_b_periodic(b, rho, dt, u, v, dx, dy)
        p = poisson_periodic(p, dx, dy, b, nit)
        
        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * dt / dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))) + F*dt

        v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
        
        # Calculate Periodic BC's

        # Periodic at x = 0

        u[1:-1, 0] = (un[1:-1, 0]-
                         un[1:-1, 0] * dt / dx *
                        (un[1:-1, 0] - un[1:-1, -1]) -
                         vn[1:-1, 0] * dt / dy *
                        (un[1:-1, 0] - un[0:-2, 0]) -
                         dt / (2 * rho * dx) * (p[1:-1, 1] - p[1:-1, -1]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
                         dt / dy**2 *
                        (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0]))) + F*dt

        v[1:-1,0] = (vn[1:-1, 0] -
                        un[1:-1, 0] * dt / dx *
                       (vn[1:-1, 0] - vn[1:-1, -1]) -
                        vn[1:-1, 0] * dt / dy *
                       (vn[1:-1, 0] - vn[0:-2, 0]) -
                        dt / (2 * rho * dy) * (p[2:, 1] - p[0:-2, -1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
                        dt / dy**2 *
                       (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))
        
        # Peridoic at x = 2

        u[1:-1, -1] = (un[1:-1, -1]-
                         un[1:-1, -1] * dt / dx *
                        (un[1:-1, -1] - un[1:-1, -2]) -
                         vn[1:-1, -1] * dt / dy *
                        (un[1:-1, -1] - un[0:-2, -1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 0] - p[1:-1, -2]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 0] - 2 * un[1:-1, -1] + un[1:-1, -2]) +
                         dt / dy**2 *
                        (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1]))) + F*dt

        v[1:-1,-1] = (vn[1:-1, -1] -
                        un[1:-1, -1] * dt / dx *
                       (vn[1:-1, -1] - vn[1:-1, -2]) -
                        vn[1:-1, -1] * dt / dy *
                       (vn[1:-1, -1] - vn[0:-2, -1]) -
                        dt / (2 * rho * dy) * (p[2:, 0] - p[0:-2, -2]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
                        dt / dy**2 *
                       (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

        u[0, :]  = 0
        u[-1, :] = 0  
        v[0, :]  = 0
        v[-1, :] = 0

        udiff = (np.sum(u)-np.sum(un))/np.sum(u)
        stepcount += 1
        
        
        
    return u, v, p

if __name__ == '__main__':
    print('zaza')